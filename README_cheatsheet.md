# photo.py Cheat Sheet (10 Zeilen)
1. Start: python photo.py /pfad/zum/ordner
2. Bestimmte Dateien: python photo.py /pfad --indices 0 1
3. Bildplot speichern: python photo.py /pfad --save frames.png
4. Adaptive Apertur je Stern: python photo.py /pfad --photometry-method aperture --aperture-per-star
5. Nahe Sterne besser trennen: python photo.py /pfad --fwhm 3.0 --detect-fwhm-scale 0.75
6. PSF (FWHM automatisch pro Bild, schlechte Fits entfernt): python photo.py /pfad --photometry-method psf [--psf-fwhm 4.0]
7. CMD speichern: python photo.py /pfad --save-cm cmd.png
8. Methodenvergleich: python photo.py /pfad --compare-cmd-methods --save-cmd-compare cmp.png
9. Sternprofil+PSF Plot: python photo.py /pfad --plot-stellar-profiles --save-stellar-profiles profile.png
10. Kalibrieren + Isochrone: python photo.py /pfad --calibrate --fit-isochrone --ebv 0.05
