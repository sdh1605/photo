# photo.py Cheat Sheet (10 Zeilen)
1. Start: python photo.py /pfad/zum/ordner
2. Bestimmte Dateien: python photo.py /pfad --indices 0 1
3. Bildplot speichern: python photo.py /pfad --save frames.png
4. Apertur: python photo.py /pfad --photometry-method aperture
5. PSF: python photo.py /pfad --photometry-method psf
6. CMD speichern: python photo.py /pfad --save-cm cmd.png
7. Methodenvergleich: python photo.py /pfad --compare-cmd-methods --save-cmd-compare cmp.png
8. Header ignorieren: python photo.py /pfad --filters B V --force-filters
9. Kalibrieren + Isochrone: python photo.py /pfad --calibrate --fit-isochrone --ebv 0.05
10. Keine Matches/Farbausreisser: --maxsep anheben, --fwhm/--threshold anpassen, --outlier-sigma erhoehen
