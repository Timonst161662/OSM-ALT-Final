# A_Star_and_ALT Routing

## Übersicht
C++17 Routing-Server mit Dijkstra, A*, ALT (uni, bi, OR).  
HTTP-Server via **Crow** auf Port `18080`.  
Antwortet auf `GET /route?source=<id>&target=<id>`.

##Testdaten
0831_1kRandomQueries_Log.txt enthält den gesamten Log der 1000 random Queries
0831_1kRandomQueries_summary.txt enthält die Zusammmenfassung der 1000 random Queries

## Voraussetzungen
- Visual Studio 2019/2022 (x64, Release, C++17)
- Abhängigkeiten (liegen im Projekt):
  - [Crow](https://github.com/CrowCpp/Crow) (header-only)
  - [nlohmann/json](https://github.com/nlohmann/json) (header-only)
- Wenn die Landmark Tabellen (tables_best_K220) als zip oder mehrere zips gedownloadet wird müssen sie entpackt in einem Ordner "tabels_best_k220" in A-Star_and_ALT abgelegt werden
- Graph und Landmark Tabellen sind zu groß für github.

## Projektstruktur
```
A-Star_and_ALT/
├─ ALT.cpp / ALT.h
├─ AStar.cpp / AStar.h
├─ Dijkstra.cpp / Dijkstra.h
├─ GraphUtils.cpp / GraphUtils.h
├─ Geo.cpp / Geo.h
├─ routing.cpp          # enthält main(...)
├─ externe includes/
│   ├─ crow_x64-windows/include/crow/app.h
│   └─ nlohmann-json_x64-windows/include/nlohmann/json.hpp
├─ tables_best_K220/    # Landmark-Tabellen (L_fwd_*.tbl, …)
├─ graph_output_fixed.fmi
└─ static/              # Frontend-Dateien (index.html, js, css)
```

##Windows
1. alles downloaden und entpacken wenn nötig
2. bash: cl /std:c++17 /EHsc /O2 ^
  /I"externe includes\crow_x64-windows\include" ^
  /I"externe includes\nlohmann-json_x64-windows\include" ^
  /I"externe includes\asio_x64-windows\include" ^
  routing.cpp ALT.cpp AStar.cpp Dijkstra.cpp GraphUtils.cpp Geo.cpp ^
  Ws2_32.lib
(3. für MinGW g++ -std=c++17 -O2 ^
  -I"externe includes\crow_x64-windows\include" ^
  -I"externe includes\nlohmann-json_x64-windows\include" ^
  -I"externe includes\asio_x64-windows\include" ^
  routing.cpp ALT.cpp AStar.cpp Dijkstra.cpp GraphUtils.cpp Geo.cpp ^
  -lws2_32 -o router.exe) --> start über router.exe


#Mac
1. alles downloaden und entpacken wenn nötig
2. bash:clang++ -std=c++17 -O2 \
  -I"externe_includes/crow_x64-windows/include" \
  -I"externe_includes/nlohmann-json_x64-windows/include" \
  -I"externe_includes/asio_x64-windows/include" \
  routing.cpp ALT.cpp AStar.cpp Dijkstra.cpp GraphUtils.cpp Geo.cpp \
  -o router
3. starte: ./router


## Start
1. Stelle sicher, dass `graph_output_fixed.fmi` und `tables_best_K220/` im Projektordner liegen.  
2. Programm starten:
   ```
   x64\Release\A_Star_and_ALT.exe
   ```
3. Browser oder `curl`:
   ```
   http://127.0.0.1:18080/route?source=12345&target=67890
   ```

## Hinweise
- Entscheidung uni/bidirectional (OR) in `ALT.cpp`:
  ```cpp
  double OR_SWITCH_METERS = 674087.637;
  ```
- `static/` enthält HTML/JS/CSS fürs Frontend.  
- Alle Pfade sind relativ → Projekt lässt sich ohne Anpassung verschieben.
- routing.cpp enthöt zwei main funktionen, die Erste führt alle 5 Algorithmen aus, die Zweite nur ALt_UniOrBi. Je nach Bedarf muss die jeweilige Funktion alleine den Namen main tragen. Die andere kann bspw. main_ genannt werden.
