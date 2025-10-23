# Course Website Template (Quarto)

Ein Quarto-Gerüst für Kurs-Webseiten mit **Branding**, **Impressum** und **Auto-Deploy auf GitHub Pages**.

Ein Blick auf die [Demo Seite](https://gisma-courses.github.io/web-template/) zeigt das Default-setup.

> Das Template kann grundsätzlich auch vollständig im Browser genutzt werden. Dann müssen jedoch alle Einträge, die `scripts/configure.py` nutzt, **online** in der `site-config.yaml` erfolgen . Das commiten löst dann die automatische Einrichtung der Seite aus
>
> Empfehlenswert ist dennoch: Repo einmal lokal klonen, die `site-config.yaml` ausfüllen und anschließend `python3 scripts/configure.py` ausführen. Danach pushen – ab dann läuft der Build/Deploy automatisch und die Inhalte können bei Bedarf auch ausschließlich online weiterbearbeitet werden.

------------------------------------------------------------------------

## Quick Setup

-   **Repo anlegen**
    -   GitHub: *Use this template* → neues Repo

    -   CLI (leerem Ordner):

        ``` bash
        quarto use template ORG/REPO
        # Trust? Y   •   Create subdirectory? n
        ```
-   **Konfiguration**
    -   **Online:** `site-config.yaml` im Webeditor ausfüllen → Push.\
        *Der Workflow führt automatisch `python3 scripts/configure.py --noninteractive` aus.*

    -   **Lokal (komfortabel):**

        ``` bash
        git clone https://github.com/<USER>/<REPO>.git
        cd <REPO>
        python3 scripts/configure.py --interactive   # fragt fehlende Felder
        git commit -am "Apply site-config" && git push
        ```
-   **Build & Deploy (automatisch)**
    -   GitHub Actions rendert nach **`docs/`** und veröffentlicht via Pages.
-   **Einmalige GitHub-Einstellungen**
    -   **Settings → Pages:** *Deploy from a branch* → **Branch:** `main` • **Folder:** `/docs`
    -   **Settings → Actions → General:** *Workflow permissions* → **Read and write permissions**
-   **Tägliche Nutzung**
    -   Inhalte ändern (`*.qmd`, Bilder) → `git push` → CI baut
    -   Branding/Impressum ändern? → `site-config.yaml` anpassen\
        (lokal `python3 scripts/configure.py` **oder** nur online; Workflow setzt’s automatisch)
-   **Hinweise / Stolpersteine**
    -   Pflichtfelder in `site-config.yaml` nicht leer lassen (sonst Build-Fehler mit Hinweis)
    -   Navbar-Rechts (`portal_text`, `portal_url`) wird nur dann gesetzt/überschrieben, wenn beide Werte vorhanden sind – also portal_text und portal_url!
    -   404? Prüfe Pages-Einstellung (main/docs) & ob `docs/index.html` existiert
    -   Alles in `./template/` gelandet? Beim `quarto use` **„Create subdirectory?“ → n**

------------------------------------------------------------------------

## Voraussetzungen

-   Quarto CLI (empfohlen ≥ 1.5): `quarto --version`
-   Python 3 zum Setup: `python3 --version`
-   GitHub-Repo mit Schreibrechten

------------------------------------------------------------------------
