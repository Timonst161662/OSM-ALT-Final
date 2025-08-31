#include "Geo.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <algorithm>

namespace geo {

    // --------- Winkel/Projektion ---------
    double deg2rad(double d) { return d * PI / 180.0; }
    double rad2deg(double r) { return r * 180.0 / PI; }

    double normLonDeg(double lon) {
        double x = std::fmod(lon + 540.0, 360.0);
        if (x < 0) x += 360.0;
        return x - 180.0;
    }

    double normalizeDegrees(double deg) {
        while (deg > 180.0) deg -= 360.0;
        while (deg < -180.0) deg += 360.0;
        return deg;
    }

    double dlon_deg_wrap(double lon1, double lon2) {
        double d = lon1 - lon2;
        while (d > 180.0) d -= 360.0;
        while (d < -180.0) d += 360.0;
        return d;
    }

    double m_per_deg_lat() { return 111000.0; }
    double m_per_deg_lon_at_lat(double lat_deg) { return 111000.0 * std::cos(deg2rad(lat_deg)); }

   
    // --------- Distanzen ---------
    double haversine(double lat1, double lon1, double lat2, double lon2) {
        const double dLat = deg2rad(lat2 - lat1);
        const double dLon = deg2rad(dlon_deg_wrap(lon2, lon1));
        const double s1 = std::sin(dLat * 0.5);
        const double s2 = std::sin(dLon * 0.5);
        const double a = s1 * s1 +
            std::cos(deg2rad(lat1)) * std::cos(deg2rad(lat2)) * s2 * s2;
        return 2.0 * EARTH_RADIUS_M * std::atan2(std::sqrt(a), std::sqrt(std::max(0.0, 1.0 - a)));
    }

    // --------- AABB / Polygone / Küste ---------
    bool aabbOverlap(double aMinLon, double aMinLat, double aMaxLon, double aMaxLat,
        double bMinLon, double bMinLat, double bMaxLon, double bMaxLat)
    {
        // Lat: normaler Intervall-Test
        if (aMaxLat < bMinLat || bMaxLat < aMinLat) return false;

        // Lon: wrap-sicher
        if (!lonOverlap(aMinLon, aMaxLon, bMinLon, bMaxLon)) return false;

        return true;
    }


    bool pointNearBBox(const BBox& b, double lat, double lon, double marginDeg)
    {
        const double latMin = b.minLat - marginDeg;
        const double latMax = b.maxLat + marginDeg;
        const double lonMin = b.minLon - marginDeg;
        const double lonMax = b.maxLon + marginDeg;

        const bool latOk = (lat >= latMin && lat <= latMax);
        const bool lonOk = lonContains(lonMin, lonMax, lon);

        return latOk && lonOk;
    }

    std::vector<int> queryRTreeIndicesAABB(
        const std::vector<RTreeEntry>& rtree,
        double minLon, double minLat,
        double maxLon, double maxLat)
    {
        // Variablendeklaration:
        bool cellLonWrapSuspicious = (minLon > maxLon) || (std::abs(maxLon - minLon) > 350.0);

        std::vector<int> out;
        out.reserve(16);

        for (int i = 0; i < static_cast<int>(rtree.size()); ++i) {
            const auto& e = rtree[(size_t)i];
            if (aabbOverlap(minLon, minLat, maxLon, maxLat,
                e.bbox.minLon, e.bbox.minLat, e.bbox.maxLon, e.bbox.maxLat))
            {
                /*if (out.size() <= 3) {
                    std::cout << "[AABB:HIT] bboxLon=[" << e.bbox.minLon << "," << e.bbox.maxLon
                        << "] bboxLat=[" << e.bbox.minLat << "," << e.bbox.maxLat
                        << "] idx=" << i << "\n";
                }*/

                out.push_back(i);
            }
        }

        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        if (out.empty() && cellLonWrapSuspicious) {
            std::cout << "[AABB:WRAP?] cellLon=[" << minLon << "," << maxLon
                << "] cellLat=[" << minLat << "," << maxLat << "] -> no overlaps (suspicious)\n";
        }

        return out;
    }




    double distancePointToBBoxMeters(const Point& p, const BBox& b)
    {
        // Lon/Lat "Inside" Tests (wrap-sicher für Lon)
        const bool insideLon = lonContains(b.minLon, b.maxLon, p.x);
        const bool insideLat = (p.y >= b.minLat && p.y <= b.maxLat);
        if (insideLon && insideLat) return 0.0;

        // Lat-Klammer (normal)
        const double lat_clamp = std::clamp(p.y, b.minLat, b.maxLat);

        // Lon-Klammer (wrap-sicher): wähle die NÄCHSTE Rand-Longitude
        // Kandidaten sind die beiden Grenzen; nimm die mit minimaler |dlon|
        const double dToMin = std::abs(dlon_deg_wrap(p.x, b.minLon));
        const double dToMax = std::abs(dlon_deg_wrap(p.x, b.maxLon));
        const double lon_clamp = (dToMin <= dToMax) ? b.minLon : b.maxLon;

        if (insideLat && !insideLon) {
            // Horizontaler Abstand zur nächsten Lon-Grenze
            const double dx_deg = std::abs(dlon_deg_wrap(p.x, lon_clamp));
            return dx_deg * m_per_deg_lon_at_lat(p.y);
        }
        if (!insideLat && insideLon) {
            // Vertikaler Abstand zur nächsten Lat-Grenze
            const double dy_deg = std::abs(p.y - lat_clamp);
            return dy_deg * m_per_deg_lat();
        }

        // Ecke: Haversine zum geklammerten Eckpunkt
        return haversine(p.y, p.x, lat_clamp, lon_clamp);
    }


    double distancePointToPolygonMeters(const Point& p, const std::vector<Polygon>& multi)
    {
        // lokales equirectangular um p
        const double kx = m_per_deg_lon_at_lat(p.y);
        const double ky = m_per_deg_lat();
        auto projX = [&](double lon) { return dlon_deg_wrap(lon, p.x) * kx; };
        auto projY = [&](double lat) { return (lat - p.y) * ky; };

        auto segDist = [&](double x0, double y0, double x1, double y1, double x2, double y2) {
            const double vx = x2 - x1, vy = y2 - y1;
            const double wx = x0 - x1, wy = y0 - y1;
            const double vv = vx * vx + vy * vy;
            double t = vv > 0 ? (wx * vx + wy * vy) / vv : 0.0;
            t = std::clamp(t, 0.0, 1.0);
            const double px = x1 + t * vx;
            const double py = y1 + t * vy;
            const double dx = x0 - px;
            const double dy = y0 - py;
            return std::sqrt(dx * dx + dy * dy);
            };

        const double x0 = 0.0, y0 = 0.0; // p projiziert ist Ursprung

        double best = std::numeric_limits<double>::infinity();
        for (const auto& poly : multi) {
            if (poly.size() < 2) continue;
            for (size_t i = 0, n = poly.size(); i < n; ++i) {
                const Point& A = poly[i];
                const Point& B = poly[(i + 1) % n];
                const double x1 = projX(A.x), y1 = projY(A.y);
                const double x2 = projX(B.x), y2 = projY(B.y);
                best = std::min(best, segDist(x0, y0, x1, y1, x2, y2));
            }
        }
        if (!std::isfinite(best)) return 5e6; // fallback
        return best;
    }

    std::vector<Point> sampleCellTestPoints(double minLon, double minLat, double maxLon, double maxLat) {
        // Referenzfenster um den Zellenmittelpunkt bauen
        const double cLon = normLonDeg((minLon + maxLon) * 0.5);
        const double cLat = (minLat + maxLat) * 0.5;

        auto wrapToCenter = [&](double L) {
            // bringe L in das Fenster um cLon (kleinster Δ)
            return normLonDeg(cLon + dlon_deg_wrap(L, cLon));
            };

        return {
            Point{ wrapToCenter((minLon + maxLon) * 0.5), cLat }, // Mitte
            Point{ wrapToCenter(minLon), minLat },
            Point{ wrapToCenter(minLon), maxLat },
            Point{ wrapToCenter(maxLon), minLat },
            Point{ wrapToCenter(maxLon), maxLat }
        };
    }


    std::pair<Point, double> closestPointOnSegment(const Point& A_in, const Point& B_in, const Point& P)
    {
        // Unwrap Longitudes relativ zu P.x, damit ΔLon klein ist
        Point A{ normLonDeg(P.x + dlon_deg_wrap(A_in.x, P.x)), A_in.y };
        Point B{ normLonDeg(P.x + dlon_deg_wrap(B_in.x, P.x)), B_in.y };

        const double ax = A.x, ay = A.y;
        const double bx = B.x, by = B.y;
        const double px = P.x, py = P.y;

        const double vx = bx - ax, vy = by - ay;
        const double wx = px - ax, wy = py - ay;

        const double vv = vx * vx + vy * vy;
        double t = vv > 0.0 ? (vx * wx + vy * wy) / vv : 0.0;
        t = std::clamp(t, 0.0, 1.0);

        Point foot{ ax + t * vx, ay + t * vy };
        double d = haversine(py, px, foot.y, foot.x); // Meter
        return { foot, d };
    }


    Point advancePointMeters(const Point& S, double angRad, double meters)
    {
        const double dLatDeg = (meters * std::sin(angRad)) / 111000.0;
        const double latMid = S.y + 0.5 * dLatDeg;                 // Stabilisierung
        const double coslat = std::max(0.15, std::cos(deg2rad(latMid)));
        const double dLonDeg = (meters * std::cos(angRad)) / (111000.0 * coslat);
        return Point{ normLonDeg(S.x + dLonDeg), std::clamp(S.y + dLatDeg, -90.0, 90.0) };
    }

    

} // namespace geo

// Punkt-in-Rechteck-Test
bool Rect::contains(const geo::Point& p) const {
    return p.y >= y_min && p.y <= y_max &&
        p.x >= x_min && p.x <= x_max;
}

bool BBox::contains(const geo::Point& p) const {
    if (p.y < minLat || p.y > maxLat) return false;
    return geo::lonContains(minLon, maxLon, p.x); // wrap-sicher
}

bool BBox::intersects(const BBox& other) const {
    if (maxLat < other.minLat || minLat > other.maxLat) return false;
    return geo::lonOverlap(minLon, maxLon, other.minLon, other.maxLon);
}