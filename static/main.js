// main.js — Unendliche Tiles + route ohne Sprünge oder Lücke
/*
let map = L.map('map', {
    worldCopyJump: true,
    preferCanvas: true
}).setView([48.7758, 9.1829], 3);

L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; OSM contributors',
    noWrap: false // Tiles dürfen wrappen (Standard)
}).addTo(map);

let nodeMarkers = {};
let edgeLines = [];
let pathLayer = null;         // LayerGroup der aktuell gerenderten Kopie
let pathUnwrapped = null;     // entrollter Pfad (Array [lat, lon_cont])
let pathCopies = [];          // vorbereitete Kopien (Arrays von LatLng)
let lastCopyIndex = null;

// ---- Utilities ---------------------------------------------------------

// 1) Längengrade entrollen: aufeinanderfolgende Punkte max ±180 auseinander
function unwrapPathNodes(nodes) {
    if (!nodes || nodes.length === 0) return [];
    const out = [];
    let prevLon = nodes[0].lon;
    out.push([nodes[0].lat, prevLon]);
    for (let i = 1; i < nodes.length; i++) {
        const lat = nodes[i].lat;
        let lon = nodes[i].lon;
        while (lon - prevLon > 180) lon -= 360;
        while (lon - prevLon < -180) lon += 360;
        out.push([lat, lon]);
        prevLon = lon;
    }
    return out;
}

// 2) Erzeuge Kopien des Pfads, um ±360 verschoben
function makeWrappedCopies(latlngs, baseShiftK) {
    // Wir bauen drei Kopien: k-1, k, k+1
    const copies = [];
    for (let dk = -1; dk <= 1; dk++) {
        const k = baseShiftK + dk;
        const shifted = latlngs.map(([lat, lon]) => [lat, lon + 360 * k]);
        copies.push(shifted);
    }
    return copies; // Index 0 => k-1, 1 => k, 2 => k+1
}

// 3) Wähle die Kopie, deren mittlerer Lon dem Karten-Center am nächsten ist
function pickBestCopyIndex(copies, centerLon) {
    let bestIdx = 0;
    let bestDist = Infinity;
    for (let i = 0; i < copies.length; i++) {
        const c = copies[i];
        const midLon = c[Math.floor(c.length / 2)][1]; // grobe Mitte reicht
        const dist = Math.abs(centerLon - midLon);
        if (dist < bestDist) { bestDist = dist; bestIdx = i; }
    }
    return bestIdx;
}

// 4) Render eine Kopie als einzelne Polyline
function renderCopy(copyLatLngs) {
    if (pathLayer) { map.removeLayer(pathLayer); pathLayer = null; }
    const line = L.polyline(copyLatLngs, { color: 'red', weight: 4, noClip: true });
    pathLayer = L.layerGroup([line]).addTo(map);
    map.fitBounds(L.latLngBounds(copyLatLngs));
}

// 5) Beim Panning ggf. andere Kopie wählen (damit Linie nie am Kartenrand „abreißt“)
function maybeRerenderForCenter() {
    if (!pathUnwrapped || pathUnwrapped.length < 2) return;
    const centerLon = map.getCenter().lng;
    // Falls noch keine Kopien existieren: initialisieren
    if (pathCopies.length === 0) {
        const firstLon = pathUnwrapped[0][1];
        const k0 = Math.round((centerLon - firstLon) / 360);
        pathCopies = makeWrappedCopies(pathUnwrapped, k0);
    }
    const idx = pickBestCopyIndex(pathCopies, centerLon);
    if (idx !== lastCopyIndex) {
        lastCopyIndex = idx;
        renderCopy(pathCopies[idx]);
    }
}

// ---- Routing & Zeichnen -----------------------------------------------

async function calculateRoute() {
    const source = document.getElementById("sourceNode").value;
    const target = document.getElementById("targetNode").value;
    if (!source || !target) {
        alert("Please enter both source and target node IDs.");
        return;
    }

    const res = await fetch(`/route?source=${source}&target=${target}`);
    const data = await res.json();

    console.log("Antwort:", data);
    if (!data.found) {
        console.log("Route nicht gefunden.");
        if (pathLayer) { map.removeLayer(pathLayer); pathLayer = null; }
        pathUnwrapped = null;
        pathCopies = [];
        lastCopyIndex = null;
        return;
    }

    console.log("Pfad-Knoten:", data.path_nodes);

    // 1) Entrollen (keine Sprünge)
    pathUnwrapped = unwrapPathNodes(data.path_nodes);

    // 2) Kopien um das aktuelle Karten-Center vorbereiten
    const centerLon = map.getCenter().lng;
    const firstLon = pathUnwrapped[0][1];
    const k0 = Math.round((centerLon - firstLon) / 360);
    pathCopies = makeWrappedCopies(pathUnwrapped, k0);

    // 3) „beste“ Kopie rendern
    lastCopyIndex = pickBestCopyIndex(pathCopies, centerLon);
    renderCopy(pathCopies[lastCopyIndex]);
}

// Bei Kartenbewegung ggf. andere Kopie anzeigen
map.on('moveend', () => {
    if (!pathUnwrapped) return;
    // Prüfen, ob Center so weit sprang, dass eine andere Kopie schöner ist
    maybeRerenderForCenter();
});

// Deine bestehende Funktion (Knoten/Edges etc. laden)
loadGraph();
*/

// main.js — Unendliche Tiles + route ohne Sprünge, plus blaue Landmark-Markierungen

let map = L.map('map', {
    worldCopyJump: true,
    preferCanvas: true
}).setView([48.7758, 9.1829], 3);

L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: '&copy; OSM contributors',
    noWrap: false // Tiles dürfen wrappen (Standard)
}).addTo(map);

// ==== State =============================================================

let nodeMarkers = {};         // wird von loadGraph() befüllt (id -> Leaflet Marker)
let edgeLines = [];
let pathLayer = null;         // LayerGroup mit aktueller Route
let pathUnwrapped = null;     // entrollter Pfad (Array [lat, lon_cont])
let pathCopies = [];          // vorbereitete Kopien (Arrays von [lat, lon])
let lastCopyIndex = null;

let landmarkLayer = null;     // LayerGroup mit Landmark-Markern (blau)

let currentWindowMidLon = null;   // Mitte (Lon) des aktuell gerenderten Routen-Fensters
let lastLmCoordsRaw = [];         // Roh-Landmark-Koordinaten (ohne Shift)


// ==== Utilities: Route ohne Sprünge / mit Wrapping ======================

function unwrapPathNodes(nodes) {
    if (!nodes || nodes.length === 0) return [];
    const out = [];
    let prevLon = nodes[0].lon;
    out.push([nodes[0].lat, prevLon]);
    for (let i = 1; i < nodes.length; i++) {
        const lat = nodes[i].lat;
        let lon = nodes[i].lon;
        while (lon - prevLon > 180) lon -= 360;
        while (lon - prevLon < -180) lon += 360;
        out.push([lat, lon]);
        prevLon = lon;
    }
    return out;
}

function shiftLonsToWindow(latlngs, midLon) {
    if (!latlngs) return [];
    return latlngs.map(([lat, lon]) => [lat, lon + 360 * Math.round((midLon - lon) / 360)]);
}

function windowMidFromCopy(copyLatLngs) {
    if (!copyLatLngs || copyLatLngs.length === 0) return null;
    return copyLatLngs[Math.floor(copyLatLngs.length / 2)][1];
}

function shiftLonsToWindow(latlngs, centerLon) {
    return latlngs.map(([lat, lon]) => {
        // bringe lon in die Nähe des Karten-Centers
        let shiftedLon = lon + 360 * Math.round((centerLon - lon) / 360);
        return [lat, shiftedLon];
    });
}

function makeWrappedCopies(latlngs, baseShiftK) {
    const copies = [];
    for (let dk = -1; dk <= 1; dk++) {
        const k = baseShiftK + dk;
        const shifted = latlngs.map(([lat, lon]) => [lat, lon + 360 * k]);
        copies.push(shifted);
    }
    return copies; // [k-1, k, k+1]
}

function pickBestCopyIndex(copies, centerLon) {
    let bestIdx = 0, bestDist = Infinity;
    for (let i = 0; i < copies.length; i++) {
        const c = copies[i];
        const midLon = c[Math.floor(c.length / 2)][1];
        const dist = Math.abs(centerLon - midLon);
        if (dist < bestDist) { bestDist = dist; bestIdx = i; }
    }
    return bestIdx;
}

function renderCopy(copyLatLngs) {
    if (pathLayer) { map.removeLayer(pathLayer); pathLayer = null; }
    const line = L.polyline(copyLatLngs, { color: 'red', weight: 4, noClip: true });
    pathLayer = L.layerGroup([line]).addTo(map);
    map.fitBounds(L.latLngBounds(copyLatLngs));
}

function maybeRerenderForCenter() {
    if (!pathUnwrapped || pathUnwrapped.length < 2) return;
    const centerLon = map.getCenter().lng;

    if (pathCopies.length === 0) {
        const firstLon = pathUnwrapped[0][1];
        const k0 = Math.round((centerLon - firstLon) / 360);
        pathCopies = makeWrappedCopies(pathUnwrapped, k0);
    }

    const idx = pickBestCopyIndex(pathCopies, centerLon);
    if (idx !== lastCopyIndex) {
        lastCopyIndex = idx;
        const chosen = pathCopies[idx];
        renderCopy(chosen);

        // Fenster-Mitte aktualisieren
        currentWindowMidLon = windowMidFromCopy(chosen);

        // Landmarks im selben Fenster neu zeichnen
        if (lastLmCoordsRaw && lastLmCoordsRaw.length > 0) {
            const lmShifted = shiftLonsToWindow(lastLmCoordsRaw, currentWindowMidLon);
            renderLandmarksBlue(lmShifted);
        }
    }
}


// ==== Utilities: Landmark-Markierungen ==================================

// Koordinaten einer Landmark aus Backend-Objekt oder ID bestimmen
function latlonFromLandmarkEntry(entry) {
    // Fall A: Objekt mit lat/lon vorhanden
    if (entry && typeof entry === 'object' && 'lat' in entry && 'lon' in entry) {
        return [entry.lat, entry.lon];
    }
    // Fall B: nur ID => aus nodeMarkers holen
    const id = typeof entry === 'object' && 'id' in entry ? entry.id : entry;
    const mk = nodeMarkers[id];
    if (mk && mk.getLatLng) {
        const ll = mk.getLatLng();
        return [ll.lat, ll.lng];
    }
    return null;
}

function flattenUsedLandmarks(data) {
    // Unterstützt: used_landmarks ODER forward/backward
    let raw = [];
    if (Array.isArray(data.used_landmarks)) {
        raw = data.used_landmarks;
    } else {
        if (Array.isArray(data.used_landmarks_forward)) raw = raw.concat(data.used_landmarks_forward);
        if (Array.isArray(data.used_landmarks_backward)) raw = raw.concat(data.used_landmarks_backward);
    }
    // In Lat/Lon umwandeln + deduplizieren
    const coords = [];
    const seen = new Set();
    for (const e of raw) {
        const ll = latlonFromLandmarkEntry(e);
        if (!ll) continue;
        const key = ll[0].toFixed(6) + ',' + ll[1].toFixed(6);
        if (!seen.has(key)) { seen.add(key); coords.push(ll); }
    }
    return coords;
}

function renderLandmarksBlue(latlngs) {
    if (landmarkLayer) { map.removeLayer(landmarkLayer); landmarkLayer = null; }
    if (!latlngs || latlngs.length === 0) return;
    const markers = latlngs.map(([lat, lon]) =>
        L.circleMarker([lat, lon], {
            radius: 6,
            color: '#1e90ff',
            fillColor: '#1e90ff',
            fillOpacity: 0.95,
            weight: 2
        })
    );
    landmarkLayer = L.layerGroup(markers).addTo(map);
}

// ==== Routing ===========================================================

async function calculateRoute() {
    const source = document.getElementById("sourceNode").value;
    const target = document.getElementById("targetNode").value;
    if (!source || !target) {
        alert("Please enter both source and target node IDs.");
        return;
    }

    const res = await fetch(`/route?source=${source}&target=${target}`);
    const data = await res.json();

    console.log("Antwort:", data);
    if (!data.found) {
        console.log("Route nicht gefunden.");
        if (pathLayer) { map.removeLayer(pathLayer); pathLayer = null; }
        if (landmarkLayer) { map.removeLayer(landmarkLayer); landmarkLayer = null; }
        pathUnwrapped = null;
        pathCopies = [];
        lastCopyIndex = null;
        currentWindowMidLon = null;
        lastLmCoordsRaw = [];
        return;
    }

    console.log("Pfad-Knoten:", data.path_nodes);

    // Route entrollen -> Kopien -> beste rendern
    pathUnwrapped = unwrapPathNodes(data.path_nodes);
    const centerLon = map.getCenter().lng;   // <-- nur EINMAL deklarieren
    const firstLon = pathUnwrapped[0][1];
    const k0 = Math.round((centerLon - firstLon) / 360);
    pathCopies = makeWrappedCopies(pathUnwrapped, k0);

    lastCopyIndex = pickBestCopyIndex(pathCopies, centerLon);
    const chosen = pathCopies[lastCopyIndex];
    renderCopy(chosen);

    // Mitte des aktiven Fensters merken
    currentWindowMidLon = windowMidFromCopy(chosen);

    // === Landmarks (blau) ===
    lastLmCoordsRaw = flattenUsedLandmarks(data);                    // Roh speichern
    const lmShifted = shiftLonsToWindow(lastLmCoordsRaw, currentWindowMidLon);
    renderLandmarksBlue(lmShifted);
}


// Bei Kartenbewegung ggf. andere Kopie anzeigen
map.on('moveend', () => {
    if (!pathUnwrapped) return;
    maybeRerenderForCenter();
});la

// Vorhandenes Laden deiner Graph-Daten (füllt nodeMarkers etc.)
loadGraph();
