/*
 * pluto-hi-monitor - Real-time Hydrogen Line (1420.405 MHz) Terminal Monitor
 * Copyright (C) 2026 Thierry VIEIL
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include <csignal>
#include <chrono>

#include <ad9361.h>
#include <iio.h>
#include <fftw3.h>

// Intégration de la bibliothèque C sigutils
extern "C" {
    #include <sigutils/sigutils.h>
    #include <sigutils/version.h>
    #include <sigutils/detect.h>
}

static volatile sig_atomic_t execution_active = 1;
void interrupteur_signal(int)
{
    execution_active = 0;
}

std::string detecter_uri_usb()
{
    std::string uri_trouvee = "";
    struct iio_scan_context *scan_ctx = iio_create_scan_context(nullptr, 0);
    if (!scan_ctx) return uri_trouvee;

    struct iio_context_info **info_list = nullptr;
    ssize_t nb = iio_scan_context_get_info_list(scan_ctx, &info_list);

    for (ssize_t i = 0; i < nb; ++i)
    {
        const char *uri = iio_context_info_get_uri(info_list[i]);
        if (uri && std::strncmp(uri, "usb:", 4) == 0)
        {
            uri_trouvee = uri;
            break;
        }
    }
    iio_context_info_list_free(info_list);
    iio_scan_context_destroy(scan_ctx);
    return uri_trouvee;
}

int main(int argc, char *argv[])
{
    std::signal(SIGINT, interrupteur_signal);

    std::string uri = (argc > 1) ? argv[1] : detecter_uri_usb();
    if (uri.empty())
    {
        std::cout << "Aucun équipement SDR USB détecté." << std::endl;
        return 0;
    }

    struct iio_context *ctx = iio_create_context_from_uri(uri.c_str());
    if (!ctx) return 1;

    const long long freq_cible_hz = 1420405000LL;
    const long long offset_lo_hz  = 1000000LL;
    const long long freq_lo_hz    = freq_cible_hz - offset_lo_hz;
    const unsigned long sample_rate = 2500000;

    struct iio_device *phy = iio_context_find_device(ctx, "ad9361-phy");
    if (!phy) { iio_context_destroy(ctx); return 1; }

    struct iio_channel *lo_rx = iio_device_find_channel(phy, "altvoltage0", true);
    if (lo_rx) iio_channel_attr_write_longlong(lo_rx, "frequency", freq_lo_hz);
    ad9361_set_bb_rate(phy, sample_rate);

    struct iio_device *rx = iio_context_find_device(ctx, "cf-ad9361-lpc");
    if (!rx) { iio_context_destroy(ctx); return 1; }

    struct iio_channel *rx_i = iio_device_find_channel(rx, "voltage0", false);
    struct iio_channel *rx_q = iio_device_find_channel(rx, "voltage1", false);
    if (!rx_i || !rx_q) { iio_context_destroy(ctx); return 1; }

    iio_channel_enable(rx_i);
    iio_channel_enable(rx_q);

    const int N = 4096;
    struct iio_buffer *buf = iio_device_create_buffer(rx, N, false);
    if (!buf) { iio_context_destroy(ctx); return 1; }

    fftw_complex *in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    std::vector<double> fenetre(N);
    for (int i = 0; i < N; ++i)
    {
        fenetre[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (N - 1)));
    }

    const int NB_LIGNES = 21;
    const int CANAL_HI_INDEX = NB_LIGNES / 2;
    const double bande_affichee_hz = 210000.0;
    const double pas_canal_hz = bande_affichee_hz / NB_LIGNES;
    const double f_debut_rf = freq_cible_hz - (bande_affichee_hz / 2.0);
    const double bin_width = (double)sample_rate / N;
    const int LARGEUR_MAX_BARRE = 38;

    struct CanalInfo {
        int debut;
        int fin;
        int count;
    };
    std::vector<CanalInfo> canaux(NB_LIGNES);
    for (int c = 0; c < NB_LIGNES; ++c)
    {
        double f_canal_rf = f_debut_rf + (c + 0.5) * pas_canal_hz;
        double f_canal_bb = f_canal_rf - freq_lo_hz;
        int bin_centre = (int)std::round(f_canal_bb / bin_width);
        int demi_largeur = (int)std::max(1.0, std::round((pas_canal_hz / 2.0) / bin_width));
        canaux[c].debut = std::max(1, bin_centre - demi_largeur);
        canaux[c].fin   = std::min(N / 2 - 1, bin_centre + demi_largeur);
        canaux[c].count = std::max(1, canaux[c].fin - canaux[c].debut + 1);
    }

    // Initialisation du détecteur sigutils (fenêtre de 16 tranches, seuil +3.0 dB)
    su_peak_detector_t detecteur_hi;
    su_peak_detector_init(&detecteur_hi, 16, 3.0);

    std::vector<double> accum_lin(NB_LIGNES, 0.0);
    std::vector<double> spectre_integre_db(NB_LIGNES, 0.0);
    std::vector<double> pics_db(NB_LIGNES, -1e9);

    const double DUREE_INTEGRATION_SEC = 15.0;
    const double DECAY_DB_PAR_FRAME = 0.08;

    const std::string RESET       = "\033[0m";
    const std::string BOLD        = "\033[1m";
    const std::string VERT        = "\033[32m";
    const std::string JAUNE       = "\033[33m";
    const std::string ROUGE       = "\033[31m";
    const std::string MAGENTA     = "\033[1;35m";
    const std::string CYAN_GRAS   = "\033[1;36m";
    const std::string BLANC_GRAS  = "\033[1;37m";
    const std::string ALERTE_FOND = "\033[41;1;37m";
    const std::string EFFACE_FIN  = "\033[K";

    std::cout << "\033[2J\033[?25l" << std::flush;

    auto t_debut_session = std::chrono::steady_clock::now();
    auto t_dernier_affichage = std::chrono::steady_clock::now();
    auto t_derniere_detection = std::chrono::steady_clock::time_point::min();

    unsigned long long nb_fft_cycle = 0;
    unsigned long long nb_fft_total = 0;
    unsigned long long total_detections = 0;

    while (execution_active)
    {
        ssize_t nbytes = iio_buffer_refill(buf);
        if (nbytes <= 0) break;

        int16_t *samples = (int16_t *)iio_buffer_start(buf);

        for (int i = 0; i < N; ++i)
        {
            in[i][0] = (double)samples[2 * i] * fenetre[i];
            in[i][1] = (double)samples[2 * i + 1] * fenetre[i];
        }

        fftw_execute(plan);
        nb_fft_cycle++;
        nb_fft_total++;

        for (int c = 0; c < NB_LIGNES; ++c)
        {
            double sum_canal = 0.0;
            for (int b = canaux[c].debut; b <= canaux[c].fin; ++b)
            {
                sum_canal += out[b][0] * out[b][0] + out[b][1] * out[b][1];
            }
            accum_lin[c] += (sum_canal / canaux[c].count);
        }

        auto t_courant = std::chrono::steady_clock::now();
        std::chrono::duration<double> ecoulee_affiche = t_courant - t_dernier_affichage;

        if (ecoulee_affiche.count() >= 0.20)
        {
            std::chrono::duration<double> ecoulee_cycle = t_courant - t_debut_session;
            double progression = std::min(1.0, ecoulee_cycle.count() / DUREE_INTEGRATION_SEC);

            double min_db = 1e9, max_db = -1e9;
            for (int c = 0; c < NB_LIGNES; ++c)
            {
                double moy_lin = accum_lin[c] / std::max(1ULL, nb_fft_cycle);
                double db = (moy_lin > 1e-12) ? (10.0 * std::log10(moy_lin)) : 0.0;
                spectre_integre_db[c] = db;

                if (db > pics_db[c])
                {
                    pics_db[c] = db;
                }
                else
                {
                    pics_db[c] -= DECAY_DB_PAR_FRAME;
                }

                if (db < min_db) min_db = db;
                if (pics_db[c] > max_db) max_db = pics_db[c];
            }

            // Détection sigutils sur le canal HI (float réel)
            double moy_hi_lin = accum_lin[CANAL_HI_INDEX] / std::max(1ULL, nb_fft_cycle);
            float ech_hi = (float)std::sqrt(moy_hi_lin);

            if (su_peak_detector_feed(&detecteur_hi, ech_hi))
            {
                total_detections++;
                t_derniere_detection = t_courant;
            }

            std::chrono::duration<double> duree_depuis_det = t_courant - t_derniere_detection;
            bool alerte_active = (duree_depuis_det.count() < 2.5);

            double delta_db = max_db - min_db;
            if (delta_db < 3.0)
            {
                max_db = min_db + 3.0;
                delta_db = 3.0;
            }

            std::cout << "\033[H";
            std::cout << BOLD << BLANC_GRAS 
                      << "====== SPECTRE INTÉGRÉ RADIOASTRONOMIE (15s) [Ctrl+C] ======" 
                      << RESET << EFFACE_FIN << "\n";
            std::cout << "LO : " << (freq_lo_hz / 1e6) << " MHz | Cible : " << (freq_cible_hz / 1e6) 
                      << " MHz (Offset : +" << (offset_lo_hz / 1e3) << " kHz)" << EFFACE_FIN << "\n";
            
            int len_prog = (int)(progression * 20);
            std::string jauge = std::string(len_prog, '=') + std::string(20 - len_prog, ' ');
            std::cout << "Cycle : [" << CYAN_GRAS << jauge << RESET << "] " 
                      << std::fixed << std::setprecision(1) << ecoulee_cycle.count() << "s / " 
                      << (int)DUREE_INTEGRATION_SEC << "s | FFT : " << nb_fft_total << EFFACE_FIN << "\n";

            std::cout << "Échelle : min=" << min_db << " dB | max=" << max_db << " dB (span " 
                      << std::setprecision(1) << delta_db << " dB) | sigutils : ";
            if (alerte_active)
            {
                std::cout << ALERTE_FOND << " !! SIGNAL HI DÉTECTÉ !! " << RESET 
                          << " (Total: " << total_detections << ")";
            }
            else
            {
                std::cout << VERT << "Veille active" << RESET 
                          << " (Total: " << total_detections << ")";
            }
            std::cout << EFFACE_FIN << "\n\n";

            for (int c = 0; c < NB_LIGNES; ++c)
            {
                double freq_mhz = (f_debut_rf + (c + 0.5) * pas_canal_hz) / 1e6;
                double val_db = spectre_integre_db[c];
                double pic_db = pics_db[c];

                double ratio_val = std::clamp((val_db - min_db) / delta_db, 0.0, 1.0);
                int len_barre = (int)(ratio_val * LARGEUR_MAX_BARRE);

                double ratio_pic = std::clamp((pic_db - min_db) / delta_db, 0.0, 1.0);
                int pos_pic = (int)(ratio_pic * (LARGEUR_MAX_BARRE - 1));

                std::string coul_barre = VERT;
                if (ratio_val >= 0.75)      coul_barre = ROUGE;
                else if (ratio_val >= 0.40) coul_barre = JAUNE;

                std::string ligne_visu = coul_barre;
                for (int i = 0; i < len_barre; ++i)
                {
                    ligne_visu += "█";
                }
                ligne_visu += RESET;

                if (pos_pic > len_barre)
                {
                    int espaces_avant = pos_pic - len_barre;
                    int espaces_apres = LARGEUR_MAX_BARRE - pos_pic - 1;
                    if (espaces_apres < 0) espaces_apres = 0;

                    ligne_visu += std::string(espaces_avant, ' ') 
                               + MAGENTA + "│" + RESET 
                               + std::string(espaces_apres, ' ');
                }
                else
                {
                    int espaces_restants = LARGEUR_MAX_BARRE - len_barre;
                    if (espaces_restants < 0) espaces_restants = 0;
                    ligne_visu += std::string(espaces_restants, ' ');
                }

                bool est_canal_hi = (c == CANAL_HI_INDEX);
                std::string tag_hi = est_canal_hi ? (CYAN_GRAS + " [HI]" + RESET) : "     ";
                std::string coul_freq = est_canal_hi ? CYAN_GRAS : "";

                std::cout << coul_freq << std::fixed << std::setprecision(3) << freq_mhz << " MHz" << RESET
                          << tag_hi << " [" << ligne_visu << "] "
                          << coul_barre << std::right << std::setprecision(1) << std::setw(5) << val_db << " dB" 
                          << RESET << EFFACE_FIN << "\n";
            }

            std::cout << BOLD << BLANC_GRAS 
                      << "============================================================" 
                      << RESET << EFFACE_FIN << "\n" << std::flush;

            t_dernier_affichage = t_courant;

            if (ecoulee_cycle.count() >= DUREE_INTEGRATION_SEC)
            {
                std::fill(accum_lin.begin(), accum_lin.end(), 0.0);
                nb_fft_cycle = 0;
                t_debut_session = std::chrono::steady_clock::now();
            }
        }
    }

    std::cout << "\033[?25h\033[0m\nArrêt du traitement et libération mémoire..." << std::endl;

    su_peak_detector_finalize(&detecteur_hi);
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    iio_buffer_destroy(buf);
    iio_context_destroy(ctx);

    return 0;
}
