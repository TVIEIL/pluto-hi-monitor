# Pluto HI Monitor

Moniteur spectral temps réel en mode console pour l'observation de la raie d'hydrogène neutre (HI à 1420,405 MHz) avec l'ADALM-PLUTO (AD9361).

<img width="1089" height="658" alt="Capture d’écran du 2026-10-25 19-30-33" src="https://github.com/user-attachments/assets/a1141201-fd9e-44c1-b9fe-a1bfd5beeb6e" />


## Fonctionnalités

- **Interface terminal ANSI :** 21 canaux spectralement intégrés sans scintillement.
- **Accélération FFT :** Calculs via FFTW3 avec fenêtrage de Hanning.
- **Détection adaptative :** Surveillance temps réel du canal HI via `su_peak_detector_t` de `sigutils`.
- **Zéro dépendance lourde :** Aucune interface graphique requise (compatible SSH et nœuds distants).

## Dépendances

- `libiio` et `libad9361-iio`
- `libfftw3-dev`
- `sigutils` (BatchDrake)

### Compilation de sigutils

```bash
git clone  https://github.com/BatchDrake/sigutils.git
cd sigutils
git submodule update --init --recursive
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/.local ..
make -j$(nproc)
make install
```

## Compilation du projet

```bash
export PKG_CONFIG_PATH=$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH
mkdir build && cd build
cmake ..
make -j$(nproc)
./spectre_integration
```
