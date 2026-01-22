; *** ChallengeOne Part 1***
; Programming Code ChallengeOne Part 1 Chris Reudenbach
;Agenten bewegen sich in einer Hügellandschaft und versuchen die Gipfel zu erreichen
; *** End of NetLogo 4.1beta3 Code Example Notice ***
turtles-own [peak? ]
patches-own [elevation]     ; jeder patch kennt seine Höhe
globals     [maxelevation peak avelev turtle-data tc]  ; für alle Routinen soll die maximale Höhe bekannt sein

;; to setup generiert die Hügelwelt und bevölkert sie mit Turtles
;; mittels Schieberegler wird eine definierbare Anzahl Hügelkuppen erzeugt
;; den patches werden Höhenwerte zugewiesen und nach Höhe von Grün über Braun farblich skaliert
;; mittels Schieberegler wird eine definierbare Anzahl Turtles erzeugt
to setup

  ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  ca ; löschen aller Variableninhalte
  reset-ticks
  ;; ziehe aus allen patches zufällig die im Slider eingestellte Anzahl von Patches
  ;; und weise ihnen zufällige Höhen zwischen 90 und 100 zu
  ask n-of anzahlhuegel patches [ set elevation (random-float 1000 + 900)]

  ;; Verteile von diesen patches aus den Höhenwert in die Nachbarschaft
  repeat (max-pxcor) [ diffuse elevation 1]

  ;; identifiziere die höchste Höhe
  set maxelevation max [elevation] of patches

  ;;skaliere die Farbe Braun zwischen dem Wert 0.5*maxelevation und maxelevation
  ;;skaliere die Farbe Grün zwischen dem Wert 0 und 0.5*maxelevation
  ask patches [ ifelse (elevation > maxelevation / 2) [set pcolor scale-color brown elevation  (maxelevation + 1) (maxelevation * .5) ]
                                            [set pcolor scale-color green elevation  -1 (maxelevation * 0.5)  ]]

  ;; Erzeuge die im Slider eingestellte Anzahl turtles auf zufälligen patches und färbe sie blau
  ;; hierzu dient der Befehl sprout
  ask n-of num-turtles patches [sprout 1 [set color blue
                                          set peak? false
                                          if pen-down? [pen-down]
                                          set heading 0
                               ]         ]

  ;; zuweisen des Grafikfensters
  setup-plot
  ;; exportieren der Turtlepositionen
  export
  ;;  schliessen der Exportdatei
  file-close
  ;; setzt die erste Farbe der importierten Turtles auf grau (=5)
  set tc 5


end
;;----------------------End setup ---

;; Endlosroutine (siehe Button) läßt die turtles nach dem höchsten Punkt suchen
;; go-turtles ist die Referenzroutine
;; go-observer8 und go-observer4 zeigen den Aufruf vom Observer für eine 8er NAchbarschaft bzw. eine 4erNAchbarschaft

to go-turtles
   ;;Falls alle Turtles oben (Achtung "oben" wird je nach Schalter unterschiedlich definiert) sind halte an
   if all? turtles [peak?]
    [ stop ]
    ;; schreibe die höhe der derzeitigen position des turtles in height
    let height [elevation] of patch-here
;;  Die Auswahl unterschiedlicher Suchstrategien wird hier über zwei ifelse abfragen umgesetzt
;;  das heißt zuerst wird der Schalter allwissen (abfragen der höhe aller patches) abgefragt ist dieser nicht gesetzt erfolgt die abfrage auf teilwissen (Befehl in-cone)
;;  ist dieser auch nicht gesetzt laufen die turtles nach dem bekannten uphill Befehl
;;  falls schalter all-wissen gesetzt ist bewegt sich das turtle zum patch mit der maximalen Höhe und setzt peak? auf wahr
  ifelse all-wissen? [ move-to max-one-of patches [elevation] set peak? true]
                              ;; Falls der schalter nicht gesetzt ist aber der schalter teil-wissen gesetzt ist
                              ;; bewegt sich das turtle zum patch mit der maximalen Höhe innerhalb de incones mit den parametern
                              ;; vision-radius vision-angle  dann wird sich zufällige gedreht falls dann die höhe gleich der maxhöhe ist wird peak? auf wahr gesetzt
                             [ ifelse teil-wissen? [move-to max-one-of patches in-cone vision-radius vision-angle [elevation] wiggle if elevation = maxelevation [ set peak? true ] ]
                                                        ;; Falls beide Schalter nicht gesetzt sind wird die Uphill Funktion verwendet und schrittweise zum Nächsten Gipfel gelaufen
                                                        [uphill elevation  if elevation = maxelevation [ set peak? true]
                                                        ]
                             ]

  ;; berechnen der mittleren Höhe Aller Turtles
  set avelev mean [elevation] of turtles
  ;; Zählen wieviele turtles schon oben sind (wird automatisch im Outputfenster das mit dieser Variablen verbunden ist ausgegeben)
  set peak count turtles with [peak?]
  ;; plotten der mittleren Höhe Aller Turtles
  plot avelev

end
;;----------------------End go-turtles ---



;; Endlosroutine (siehe Button) läßt die turtles nach dem höchsten Punkt suchen
;; ist vergleichbar mit go-turtles uphill mit dem Unterschied dass der Befehl vom OBSERVER aufgerufen wird und neighbors benutzt
to go-observer8

;; schau in der 8-er Nachbarschaft nach dem höchsten Wert und richte dich dorthin aus
 ask turtles [ face max-one-of neighbors [elevation]
               ;; gehe einen patch vorwärts
               fd 1
               ;; plotten der mittleren Höhe Aller Turtles
               plot mean [elevation] of turtles
             ]
 end

;; Endlosroutine (siehe Button) läßt die turtles nach dem höchsten Punkt suchen
;; ist identisch mit go-observer8 mit dem Unterschied dass nur 4 Nachbarn betrachtet werden

to go-observer4
;; schau in der 4-er Nachbarschaft nach dem höchsten Wert und richte dich dorthin aus
 ask turtles [ face max-one-of neighbors4 [elevation]
               fd 1
               ;; plotten der mittleren Höhe Aller Turtles
               plot mean [elevation] of turtles
             ]
end


;; weist dem Grafikfenster den Namen  und die Eigenschaften zu
to setup-plot
  set-current-plot "average elevation of turtles"
;; skaliert die y-achse auf den Wert der nach maxheight folgenden ganzen Zahl z.B.(55.8 -> 56)
  set-plot-y-range 0 ceiling maxelevation
end

;;-------------------- Tools ---------------

;; routine die die turtles zufällig links und rechts dreht um ihnen eine neue richtung zu geben
to wiggle
  ;; Rechtsdrehung zufallswinkel 0-180
  rt random 180
  lt random 180
end

To export
    ;; Falls die Datei existiert löschen wir sie weil sonst
    ;; neue Daten an die existierenden angehangen werden
    if file-exists? "firstposition.txt"[ file-delete "firstposition.txt" ]
    ;; öffnen der Datei
    file-open "firstposition.txt"
    ;; wir bitten die trtles ihre xy positionen in die Datei zu schreiben
    ask turtles [ file-write xcor file-write ycor ]
end

To import
    ;; wir lesen die daten in eine sogenannte liste ein die besteht aus tupel mit xy koordinaten
    ;; die leere liste definieren wir hier
    set turtle-data []

    ;; öffnen der Datei
    file-open "firstposition.txt"

    ;; Einelesen aller Daten aus der Datei
    while [ not file-at-end? ]
    [
      ;; file-read liest die variablen ein in diesem Fall koordinaten also Zahlen
      ;;sie werden in einer doppelliste gesichert (z.B. [[1 1] [1 2] ...])
      ;; Beim nachfolgenden Befehl wird in jedem Schritt ein Koordinatenpärchen eingelesen und als 2er Liste zusammengefasst
      ;; in turtle-data steht dann eine liste mit koordinatenpaaren
      set turtle-data sentence turtle-data (list (list file-read file-read))
    ]
    ;; fertig datei schliessen
    file-close
    ;; Löschen der alten Turtles und
    ct
    ;; generieren der neuen an der alten ausgangsposition
    ;; für jede koordinate aus turtledata lies den wert und lasse auf diesem patch turtles mit den genannten eigenschaften entstehen
    foreach turtle-data [ ?1 -> ask patch first ?1 item 1 ?1 [ sprout 1  [set color tc
                                                                    set heading 0
                                                                    set peak? false
                                                                    if pen-down? [pen-down]
                                                                 ]
                                                       ]
                         ]

    set tc tc + 10
end
@#$#@#$#@
GRAPHICS-WINDOW
300
10
749
460
-1
-1
4.8413
1
10
1
1
1
0
0
0
1
-31
31
-31
31
0
0
1
ticks
30.0

SLIDER
15
80
150
113
anzahlhuegel
anzahlhuegel
1
100
100.0
1
1
NIL
HORIZONTAL

SLIDER
15
115
150
148
num-turtles
num-turtles
0
1000
100.0
1
1
NIL
HORIZONTAL

PLOT
14
328
174
461
average elevation of turtles
Ticks
Höhe
0.0
5.0
0.0
5.0
true
false
"" ""
PENS
"default" 1.0 0 -7500403 true "" ""
"meanelevation" 2.0 0 -13345367 true "" ""

BUTTON
151
256
280
289
Alternative Suche
go-turtles
T
1
T
TURTLE
NIL
NIL
NIL
NIL
1

BUTTON
17
40
80
73
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
19
255
138
288
NIL
go-observer8
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
20
291
138
324
NIL
go-observer4
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SWITCH
15
150
150
183
pen-down?
pen-down?
0
1
-1000

SWITCH
155
80
280
113
All-Wissen?
All-Wissen?
1
1
-1000

SWITCH
156
112
281
145
Teil-Wissen?
Teil-Wissen?
1
1
-1000

SLIDER
155
180
280
213
vision-angle
vision-angle
1
360
248.0
1
1
NIL
HORIZONTAL

SLIDER
156
147
281
180
vision-radius
vision-radius
0
100
64.0
1
1
NIL
HORIZONTAL

TEXTBOX
154
226
287
282
Erweiterte Gipfelsuche\n--------------------------------
11
105.0
1

TEXTBOX
22
228
150
254
Einfache Gipfelsuche\n-------------------------------
10
27.0
1

TEXTBOX
19
10
220
34
Basis Setup für Turtles und Hügel\n--------------------------------------------------
11
0.0
1

MONITOR
174
416
291
461
GipfelTurtles
peak
17
1
11

MONITOR
175
327
291
372
Höchster Gipfel
maxelevation
1
1
11

MONITOR
174
372
291
417
Mittlere TurtleHöhe
avelev
1
1
11

BUTTON
92
41
219
74
Import Turtleposition
import
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

@#$#@#$#@
## WHAT IS IT?
ChallengeOne Part 1 
Agenten bewegen sich in einer Hügellandschaft und versuchen die Gipfel zu erreichen

## HOW IT WORKS

**setup** generiert die Hügelwelt und bevölkert sie mit Turtles, mittels Schieberegler wird eine definierbare Anzahl Hügelkuppen erzeugt und den patches werden Höhenwerte zugewiesen und nach Höhe von Grün über Braun farblich skaliert

**go-observer8**/**go-observer4** suchen höhere patches in einer 8er bzw. 4er Nachbarschaft

**Alternative Suchstrategien** Die Auswahl unterschiedlicher Suchstrategien wird wie folgt durchgeführt. 1. allwissen  ja/nein  falls nein teilwissen mit den dort vorgeneommen Einstellungen 3. falls nein  einfacher uphill Befehl

  

## HOW TO USE IT
Testen sie die unterschiedlichen einstellung aus und versucehn sie die resultierenden Muster zu erklären. 

**Import Turtleposition**  lädt die nach dem Setup automatisch gesicherte Positionen der Turtles ein. Durch mehrfaches Betätigen können unterschiedliche Farben gewählt werden.
Gedacht ist dieses Tool zur Visualisierung unterschiedlichen Raumverhaltens der Turtles gedacht

## THINGS TO NOTICE


## THINGS TO TRY

## EXTENDING THE MODEL

## NETLOGO FEATURES

## RELATED MODELS

## CREDITS AND REFERENCES
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
1
@#$#@#$#@
