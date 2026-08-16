# photo.py

CLI-Tool zur Photometrie von zwei FITS-Bildern mit anschliessendem Color-Magnitude-Diagramm (CMD).

Das Skript:
- findet zwei FITS-Dateien in einem Ordner,
- detektiert Sterne,
- fuehrt Photometrie durch (Apertur oder PSF),
- matched Sterne zwischen zwei Filtern (typisch B und V),
- erstellt CMD-Plots,
- kann kalibrieren, deredden und Isochronen anpassen,
- kann Apertur vs. PSF direkt vergleichen.

## Voraussetzungen

Benötigte Python-Pakete:
- numpy
- matplotlib
- astropy
- photutils

Beispiel-Installation:

```bash
pip install numpy matplotlib astropy photutils
```

## Schnellstart

```bash
python photo.py /pfad/zum/fits-ordner
```

Nur bestimmte Dateien (Index 0 und 3) und Bild-Plot speichern:

```bash
python photo.py /pfad/zum/fits-ordner --indices 0 3 --save frames.png
```

Mit expliziter Methode:

```bash
python photo.py /pfad/zum/fits-ordner --photometry-method aperture
python photo.py /pfad/zum/fits-ordner --photometry-method psf
```

CMD speichern:

```bash
python photo.py /pfad/zum/fits-ordner --save-cm cmd.png
```

Sternprofile (hell/mittel/dunkel) plus PSF-Fit anzeigen:

```bash
python photo.py /pfad/zum/fits-ordner --plot-stellar-profiles
```

Vergleich Apertur vs. PSF als gemeinsamer CMD-Plot:

```bash
python photo.py /pfad/zum/fits-ordner --compare-cmd-methods --save-cmd-compare cmd_compare.png
```

## Wichtige Optionen

- `--photometry-method {aperture,psf}`
  - Waehlt die Photometrie-Methode.
- `--indices, -i`
  - Welche FITS-Dateien (0-basiert) genutzt werden. Ohne Angabe: die ersten zwei.
- `--save, -s`
  - Speichert den Bildplot der FITS-Frames als PNG.
- `--aperture FLOAT`
  - Aperturradius in Pixeln. Ohne Angabe: automatische Wahl (empfohlen).
- `--aperture-per-star`
  - Passt Apertur/Annulus pro Stern individuell an.
  - Wenn `--aperture` gesetzt ist, dient dieser Wert als Basis und wird pro Stern skaliert.
- `--fwhm FLOAT`
  - FWHM fuer Quellenfindung.
- `--detect-fwhm-scale FLOAT`
  - Skaliert die FWHM nur fuer die Peak-Erkennung.
  - Kleinere Werte trennen nahe Sterne besser (z. B. 0.75), koennen aber mehr Fehlpeaks erzeugen.
- `--threshold FLOAT`
  - Sigma-Schwelle fuer Quellenfindung.
- `--maxsep FLOAT`
  - Maximaler Matching-Abstand in Pixeln.
- `--filters B V`
  - Filterzuordnung, falls keine gute Header-Info vorliegt.
- `--force-filters`
  - Erzwingt die `--filters`-Zuordnung.
- `--outlier-sigma FLOAT`
  - Sigma-Schwelle fuer Ausreisser in B-V (Standard: 3.0).
- `--cluster-radius`, `--cluster-center`, `--annulus-radii`
  - Raeumliche Selektion von Sternen.
- `--calibrate`, `--calibrate-file`
  - Interaktive Referenzstern-Kalibrierung laden/speichern.
  - Es koennen mehrere Referenzsterne per Klick ausgewaehlt werden.
  - Fuer jeden Referenzstern werden B- und V-Magnituden abgefragt.
  - Die Kalibrierung verwendet dann den Median der einzelnen Nullpunkte.
- `--fit-isochrone`, `--ebv`, `--age-range`, `--age`
  - Isochronen-Anpassung fuer Entfernung/Alter.
- `--compare-cmd-methods`, `--save-cmd-compare`
  - Direkter Methodenvergleich im CMD.
- `--plot-stellar-profiles`, `--save-stellar-profiles`
  - Plot fuer radial gemittelte Sternprofile (hell/mittel/dunkel) inklusive angepasster PSF.

## Workflow (kurz)

1. Zwei FITS-Dateien laden.
2. Sterne detektieren.
3. Photometrie pro Bild berechnen.
4. Sterne zwischen beiden Bildern robust matchen.
5. B und V zuordnen.
6. CMD erstellen (instrumentell oder kalibriert).
7. Optional: Mehrere Referenzsterne auswaehlen und kalibrieren.
8. Optional: Dereddening und Isochrone.
9. Optional: Methodenvergleich Apertur vs. PSF mit Linien und Abweichungs-Statistik.

## Mehrere Referenzsterne kalibrieren

Bei `--calibrate` kann der Nutzer nun mehrere Referenzsterne interaktiv auswaehlen:

- Bild oeffnen und beliebig viele passende Sterne per Klick markieren
- fuer jeden Klick nach B- und V-Magnitude fragen
- Messwerte aller gewaehlten Sterne in einem gemeinsamen Kalibrierungs-File speichern
- die Kalibrierung wird aus dem Median der einzelnen Nullpunkt-Loesungen bestimmt

Das Format der Kalibrierungsdatei bleibt kompatibel zum bisherigen Einzelstern-Format, erweitert aber auf mehrere Zeilen pro Datei:

```text
# Format: ref_index ref_x ref_y mag_B mag_V
12 123.45 456.78 12.345 11.234
27 200.12 410.66 12.890 11.770
```

Einzelne Referenzsterne funktionieren weiterhin genauso wie bisher.

## Unterschied: Apertur vs. PSF-Photometrie

### Apertur-Photometrie

Prinzip:
- Summiert den Sternfluss in einem Kreis um die Sternposition.
- Schaetzt den lokalen Hintergrund in einem Ring (Annulus).
- Subtrahiert Hintergrundanteil von der Apertur-Summe.

Vorteile:
- Einfach, schnell, robust.
- Gut fuer eher isolierte Sterne und moderate Felddichte.

Nachteile:
- Bei ueberlappenden Sternen wird Fluss leicht vermischt.
- Abhaengiger von Radius-Wahl und Crowding.

### PSF-Photometrie

Prinzip:
- Passt ein PSF-Modell (hier Gaussian-basiert) an Sternprofile an.
- Trennt ueberlappende Quellen besser durch Modellanpassung.
- Liefert modellbasierten Sternfluss (`flux_fit`).

Vorteile:
- Besser bei dichtem Sternfeld (Crowding).
- Bessere Trennung naher Sterne.

Nachteile:
- Rechenaufwendiger.
- Abhaengig von der Guete des PSF-Modells.
- Bei Problemen wird im Skript auf Apertur-Photometrie zurueckgefallen.

## Trennung naher Sterne

Wenn zwei helle Sterne sehr nahe beieinander liegen und als ein Peak erkannt werden:
- `--detect-fwhm-scale` reduzieren (z. B. 0.85 -> 0.75 -> 0.65).
- `--threshold` bei Bedarf leicht erhoehen, falls zu viele Falsch-Peaks entstehen.

Beispiel:

```bash
python photo.py /pfad/zum/fits-ordner --fwhm 3.0 --detect-fwhm-scale 0.75
```

## Methodenvergleich im Skript

Mit `--compare-cmd-methods` wird ein gemeinsamer CMD-Plot erzeugt:
- Apertur- und PSF-Punkte in unterschiedlichen Farben.
- Verbindungssegmente zwischen identischen Sternen.
- Kennzahlen aus den Differenzen (PSF minus Apertur):
  - RMS Delta(B-V)
  - RMS Delta(V)
  - mittlere Tendenz in Farbe (blauer/roeter)
  - mittlere Tendenz in Helligkeit (heller/dunkler)

Falls Kalibrierung aktiv ist, wird der Vergleich ebenfalls kalibriert berechnet.
Falls `E(B-V)` verfuegbar ist, erfolgt zusaetzlich dereddening im Vergleichsplot.

## Ausgabedateien

Je nach Optionen werden u.a. erzeugt:
- Bildplot der Frames
- CMD-Plot
- CMD-Methodenvergleich
- Sternprofil/PSF-Plot (`--plot-stellar-profiles`)
- `calib.txt` (Kalibrierungsdaten; auch mehrere Referenzsterne pro Datei moeglich)
- `photometry.txt` (Sternliste mit V und B-V)
- `cluster_parameters.txt` (bei Isochronen-Fit)
- `isochrone_fit.png` (bei Isochronen-Fit)

Standardmaessig werden relative Ausgabepfade im FITS-Ordner gespeichert.

## Typische Probleme

- `photutils` nicht installiert:
  - Fehlermeldung bei Import oder beim Start der Photometrie.
  - Loesung: `pip install photutils`.
- Keine Matches gefunden:
  - `--maxsep` vergroessern,
  - `--fwhm` und `--threshold` an Daten anpassen,
  - Filterzuordnung pruefen (`--filters`, `--force-filters`).
- Unerwartete Farben im CMD:
  - Kalibrierung pruefen,
  - B/V-Zuordnung pruefen,
  - outlier-Schwelle mit `--outlier-sigma` anpassen.

  Hinweis zu Defaults:
  - `--age-range` Standard ist `6.8 7.6` (log(Alter/Jahre), ca. 6-40 Myr).
