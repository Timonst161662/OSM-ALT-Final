#pragma once
#include <vector>
#include <utility>
#include <algorithm>
#include <ostream> 

namespace geo { struct Point; }

// ----- AABB-Typen -----
struct BBox {
    double minLat{}, maxLat{}, minLon{}, maxLon{};
    bool contains(const geo::Point& p) const;
    bool intersects(const BBox& other) const;
};

struct RTreeEntry {
    BBox bbox;
    std::string filePath;
};

struct Rect {
    double x_min{}, x_max{}, y_min{}, y_max{};
    bool contains(const geo::Point& p) const;
    bool intersectsCircle(const geo::Point& center, double radiusKm) const;
};


// das Projekt arbeitet in 

namespace geo {

    struct Point { double x; double y; };
    using Polygon = std::vector<Point>;

    constexpr double PI = 3.14159265358979323846;
    constexpr double EARTH_RADIUS_M = 6371000.0;

    double deg2rad(double d);
    double rad2deg(double r);
    double normLonDeg(double lon);
    double normalizeDegrees(double deg);
    double dlon_deg_wrap(double lon1, double lon2);
    double m_per_deg_lat();
    double m_per_deg_lon_at_lat(double lat_deg);

    // --- Longitude-Helfer (wrap-sicher) ---

// true, wenn [minLon, maxLon] (ggf. wrap) den Wert L abdeckt
    static inline bool lonContains(double minLon, double maxLon, double L) {
        const double lo = normLonDeg(minLon);
        const double hi = normLonDeg(maxLon);
        const double x = normLonDeg(L);
        if (lo <= hi) return (x >= lo && x <= hi);
        // Wrap-Fall: [lo,180] U [-180,hi]
        return (x >= lo || x <= hi);
    }

    // Liefert true, wenn die beiden Lon-Intervalle (jeweils ggf. wrap) sich schneiden
    static inline bool lonOverlap(double aMin, double aMax, double bMin, double bMax) {
        const double aLo = normLonDeg(aMin);
        const double aHi = normLonDeg(aMax);
        const double bLo = normLonDeg(bMin);
        const double bHi = normLonDeg(bMax);

        auto overlapNoWrap = [](double lo1, double hi1, double lo2, double hi2) {
            return !(hi1 < lo2 || hi2 < lo1);
            };

        if (aLo <= aHi && bLo <= bHi) {
            return overlapNoWrap(aLo, aHi, bLo, bHi);
        }
        if (aLo <= aHi && bLo > bHi) {
            // B wrappt: prüfe gegen beide Teilintervalle
            return overlapNoWrap(aLo, aHi, bLo, 180.0) || overlapNoWrap(aLo, aHi, -180.0, bHi);
        }
        if (aLo > aHi && bLo <= bHi) {
            // A wrappt
            return overlapNoWrap(aLo, 180.0, bLo, bHi) || overlapNoWrap(-180.0, aHi, bLo, bHi);
        }
        // Beide wrappten: dann überdeckt zumindest ein Teil
        return true; // konservativ (beide umspannen die Datumsgrenze)
    }


    // Distanzen
    double haversine(double lat1, double lon1, double lat2, double lon2);

    inline double haversine(const Point& a, const Point& b) {
        return haversine(a.y, a.x, b.y, b.x);
    }

    // ---- AABB / Polygone / Küste ----
    // Achtung: globaler Typ via :: qualifizieren!
    bool aabbOverlap(double aMinLon, double aMinLat, double aMaxLon, double aMaxLat,
        double bMinLon, double bMinLat, double bMaxLon, double bMaxLat);

    inline bool aabbOverlap(const BBox& a, const BBox& b) {
        return aabbOverlap(a.minLon, a.minLat, a.maxLon, a.maxLat,
            b.minLon, b.minLat, b.maxLon, b.maxLat);
    }

    inline bool cellOverlapsBBox(double minLon, double minLat,
        double maxLon, double maxLat,
        const BBox& b) {
        return aabbOverlap(minLon, minLat, maxLon, maxLat,
            b.minLon, b.minLat, b.maxLon, b.maxLat);
    }


    bool pointNearBBox(const ::BBox& b, double lat, double lon, double marginDeg);

    std::vector<int> queryRTreeIndicesAABB(const std::vector<::RTreeEntry>& rtree,
        double minLon, double minLat,
        double maxLon, double maxLat);

    double distancePointToBBoxMeters(const Point& p, const ::BBox& b);

    double distancePointToPolygonMeters(const Point& p, const std::vector<Polygon>& multi);

    std::vector<Point> sampleCellTestPoints(double minLon, double minLat, double maxLon, double maxLat);

    std::pair<Point, double> closestPointOnSegment(const Point& A, const Point& B, const Point& P);

    Point advancePointMeters(const Point& S, double bearingRad, double meters);

    // Longitudes ins gemeinsame Fenster um cLon ziehen (kleinstes Δ)
    static inline double wrapToCenter(double L, double cLon) {
        return normLonDeg(cLon + dlon_deg_wrap(L, cLon));
    }

    // Liefert true NUR für einen "proper" Schnittpunkt (Innenbereich beider Segmente).
  // Berührungen (Endpunkt auf Kante, Ecke getroffen, kolinear/überlappend) => false.
    static inline bool segmentsProperIntersect2D(Point A, Point B, Point C, Point D, double cLon)
    {
        // === 1. Koordinaten entwrappen (IDL/360° behandeln) ===
        A.x = wrapToCenter(A.x, cLon); B.x = wrapToCenter(B.x, cLon);
        C.x = wrapToCenter(C.x, cLon); D.x = wrapToCenter(D.x, cLon);

        // === 2. Hilfsfunktionen ===
        auto cross = [](double ax, double ay, double bx, double by) {
            return ax * by - ay * bx;
            };
        auto orient = [&](const Point& p, const Point& q, const Point& r) {
            return cross(q.x - p.x, q.y - p.y, r.x - p.x, r.y - p.y);
            };

        constexpr double eps = 1e-12;  // Toleranz

        auto oppositeSides = [&](double s1, double s2) {
            return (s1 > eps && s2 < -eps) || (s1 < -eps && s2 > eps);
            };

        // === 3. Orientierungen berechnen ===
        const double o1 = orient(A, B, C);
        const double o2 = orient(A, B, D);
        const double o3 = orient(C, D, A);
        const double o4 = orient(C, D, B);

        // === 4. Strenger Proper-Test ===
        // Nur wenn beide Segmente sich "kreuzen" und nicht nur berühren
        if (oppositeSides(o1, o2) && oppositeSides(o3, o4)) {
            return true;
        }

        // === 5. Alle Touch-/Kolinear-Fälle werden ignoriert ===
        // Früher: onSeg-Checks -> jetzt bewusst NICHT mehr,
        // weil wir nur echte Durchschnitte wollen.
        return false;
    }



    // Prüft, ob ein Punkt in einer BBox liegt
    inline bool pointInBBox(const BBox& b, double lat, double lon) {
        return lat >= b.minLat && lat <= b.maxLat &&
            geo::lonContains(b.minLon, b.maxLon, lon);
    }

    // Prüft, ob sich zwei BBoxen überschneiden
    inline bool bboxOverlap(const BBox& a, const BBox& b) {
        return !(a.maxLat < b.minLat || a.minLat > b.maxLat) &&
            geo::lonOverlap(a.minLon, a.maxLon, b.minLon, b.maxLon);
    }
} // namespace geo


