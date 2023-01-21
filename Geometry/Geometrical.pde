abstract class Geometrical {
  /*abstract Point AsPoint();
   abstract Line AsLine();
   abstract Segment AsSegment();
   abstract Range AsRange();
   abstract Circle AsCircle();
   abstract Sphere AsSphere();
   abstract Disk AsDisk();
   abstract Triangle AsTriangle();
   abstract Polygon AsPolygon();
   abstract Box AsBox();
   abstract BoundingBox AsBBox();
   abstract Rectalinear AsRectalinear();*/
}
class Line extends Geometrical
{
  Point start;
  Vector direction;
  //
  Line(Point from, Vector dir)
  {
    start = from;
    direction = dir;
  }

  Line(Point from, Normal dir)
  {
    start = from;
    direction = dir.AsVector();
  }
  Line()
  {
    start = new Point();
    direction = new Vector(1, 0);
  }

  boolean Equvialent(Line to)
  {
    if (!direction.Parallel2D(to.direction))
    {
      return false;
    }
    Vector p = start.Sub(start, to.start);
    return p.Parallel2D(to.direction);
  }

  boolean Collide2D(Line with)
  {
    if (direction.Parallel2D(with.direction))
    {
      boolean eq = Equvialent(with);
      return eq;
    }
    return true;
  }

  boolean Collide3D(Line with)
  {
    if (direction.Parallel2D(with.direction))
    {
      boolean eq = Equvialent(with);
      if (eq)
      {
        Line xz = new Line(new Point(start.x, start.z), new Vector(direction.x, direction.z));
        Line wxz = new Line(new Point(with.start.x, with.start.z), new Vector(with.direction.x, with.direction.z));
        return eq && xz.Collide2D(wxz);
      }
      return eq;
    }
    return true;
  }

  boolean Collide2D(Point p)
  {
    if (start.Collides2D(p)) { 
      return true;
    }
    Vector lp = p.Sub(p, start);
    return lp.Parallel2D(direction);
  }

  boolean Collide3D(Point p)
  {
    if (start.Collides3D(p, Epsilon)) { 
      return true;
    }
    Vector lp = p.Sub(p, start);
    return lp.Parallel3D(direction);
  }

  boolean Collide2D(Point p, double r)
  {
    if (start.Collides2D(p, r)) { 
      return true;
    }
    Vector lp = p.Sub(p, start);
    return lp.AsPoint(start.r).Parallel2D(direction, r);
  }

  boolean Collide3D(Point p, double r)
  {
    if (start.Collides3D(p, r)) { 
      return true;
    }
    Vector lp = p.Sub(p, start);
    return lp.Parallel3D(direction, r);
  }
  Line Add(Line p, Vector offset)
  {
    Point at=p.start.Add(p.start, offset).AsPoint(p.start.r);
    return new Line(at, p.direction);
  }
  Line Sub(Line p, Vector offset)
  {
    return new Line(p.start.Sub(p.start, offset).AsPoint(p.start.r), p.direction);
  }
  Point Parametric(double t)
  {
    return start.Add(start, direction.Scaled(t)).AsPoint(start.r);
  }
  Point ParametricIntersect(double t1, Line l2, double t2)
  {
    Point v1 = start.Add(start, direction.Scaled(t1)).AsPoint(start.r);
    Point v2 = l2.start.Add(l2.start, l2.direction.Scaled(t2)).AsPoint(l2.start.r);
    Point solve= new Point(v1.x*t1-v2.x*t2, v1.y*t1+v2.y*t2, v1.z*t1-v2.z*t2);
    return solve;
  }
  double Slope2D(byte axis) {
    Vector p = start.Sub(start.Unity2D(), direction.Unity2D()).AsVector();
    if (p.x != 0)
    {
      if (axis == AxisX) {
        return (p.y / p.x) * p.y;
      } else if (axis == AxisY)
      { 
        return (p.y) + p.Length2D();
      }
    } else
      if (axis == AxisY) {
        return (p.y) + p.Length2D();
      }
    return 1;
  }
};

class Range
{

  double maxim;

  double minim;

  Range(double from, double to)
  {
    minim = from; 
    maxim = to;
    Sort();
  }
  Range(Range r1, Range r2)
  {
    if (r1.minim < r2.minim) { 
      minim = r1.minim;
    } else { 
      minim = r2.minim;
    }
    if (r1.maxim > r2.maxim) { 
      maxim = r1.maxim;
    } else { 
      maxim = r2.maxim;
    }
  }
  Range()
  {
    maxim = 1;
    minim = -1;
  }
  void Sort()
  {
    if (maxim < minim)
    {
      double m = maxim;
      maxim = minim;
      minim = m;
    }
  }
  boolean Overrlapping(Range r)
  {
    Point p1 = new Point(minim, maxim);
    Point p2 = new Point(r.minim, r.maxim);
    return Overlapping2D(p1, p2, Epsilon);
  }
  boolean Overrlapping(Range r, double threshold)
  {
    Point p1 = new Point(minim, maxim);
    Point p2 = new Point(r.minim, r.maxim);
    return Overlapping2D(p1, p2, threshold);
  }
  double Clamp(double x)
  {
    if (x < minim) return minim;
    else if (x > maxim) return maxim;
    else return x;
  }
};

class Segment extends Geometrical
{
  Vector p1;
  Vector p2;

  Segment(Vector from, Vector to)
  {
    p1 = from; 
    p2 = to;
  }
  Segment()
  {
    p1 = new Vector(0, 0); 
    p2 = new Vector(1, 0);
  }
  boolean OnSide2D(Line axis)
  {
    Cartesian d1 = p1.Sub(p1.AsPoint(1), axis.start);
    Cartesian d2 = p2.Sub(p2.AsPoint(1), axis.start);
    Cartesian n = axis.direction.Rotate90Degrees2D();
    return n.DotProduct2D(n, d1) * n.DotProduct2D(n, d2) > 0;
  }
  boolean OnSide3D(Line axis)
  {
    Cartesian d1 = p1.Sub(p1.AsPoint(1), axis.start);
    Cartesian d2 = p2.Sub(p2.AsPoint(1), axis.start);
    Cartesian n = p1.CrossProduct(axis.start, axis.direction);
    return n.DotProduct(n, d1) * n.DotProduct(n, d2) > 0;
  }
  Range Project2D(Point onto)
  {
    Cartesian u = onto.Unity2D().AsVector();
    return new Range(u.DotProduct2D(u, p1), u.DotProduct2D(u, p2));
  }
  Range Project3D(Point onto)
  {
    Cartesian u = onto.Unity();
    return new Range(u.DotProduct(u, p1), u.DotProduct(u, p2));
  }
  Range Project2D(Vector onto)
  {
    Cartesian u = onto.Unity2D();
    return new Range(u.DotProduct2D(u, p1), u.DotProduct2D(u, p2));
  }
  Range Project3D(Vector onto)
  {
    Cartesian u = onto.Unity();
    return new Range(u.DotProduct(u, p1), u.DotProduct(u, p2));
  }
  boolean Collide2D(Point p)
  {
    Cartesian d = p2.Sub(p2, p1);
    Cartesian lp =p.Sub( p, p1.AsPoint(1));
    Cartesian pr = lp.Project2D(d);
    return lp == pr && pr.Length2DSquared() <= d.Length2DSquared() && pr.DotProduct2D(pr, d) >= 0;
  }
  boolean Collide3D(Point p)
  {
    Cartesian d = p2.Sub(p2, p1);
    Cartesian lp =p.Sub( p, p1.AsPoint(1));
    Cartesian pr = lp.Project(d);
    return lp == pr && pr.LengthSquared() <= d.LengthSquared() && pr.DotProduct(pr, d) >= 0;
  }
  boolean Collide2D(Point p, double r)
  {
    Cartesian d = p2.Sub(p2, p1);
    Cartesian lp =p.Sub( p, p1.AsPoint(1));
    Cartesian pr = lp.Project2D(d);
    return lp == pr && pr.Length2DSquared() <= d.Length2DSquared() && pr.DotProduct2D(pr, d) >= 0;
  }
  boolean Collide3D(Point p, double r)
  {
    Cartesian d =p2.Sub( p2, p1);
    Cartesian lp = p.Sub(p, p1.AsPoint(1));
    Cartesian pr = lp.Project(d);
    boolean ok = Math.abs(lp.x - pr.x) < r && Math.abs(lp.y - pr.y) < r && Math.abs(lp.z - pr.z) < r;
    if (ok)
    {
      double r2 = r * r;
      return pr.LengthSquared() <= (d.LengthSquared() + r2) && pr.DotProduct2D(pr, d) >= -r2;
    }
    return false;
  }
  boolean Collide3D(double thickness, Point p, double r)
  {
    Cartesian d = p2.Sub(p2, p1);
    Cartesian lp = p.Sub(p, p1.AsPoint(1));
    Cartesian pr = lp.Project(d);
    double rt = r + thickness;
    boolean ok = Math.abs(lp.x - pr.x) < rt && Math.abs(lp.y - pr.y) < rt && Math.abs(lp.z - pr.z) < rt;
    if (ok)
    {
      double r2 = rt * rt;
      return pr.LengthSquared() <= (d.LengthSquared() + r2) && pr.DotProduct2D(pr, d) >= -rt;
    }
    return false;
  }
  boolean Collide2D(Segment with)
  {
    Line ax1 = new Line();
    Line ax2 = new Line();
    ax1.start = p1.AsPoint(1);
    ax1.direction = p1.Sub(p1, p2).Unity2D().AsVector();
    if (with.OnSide2D(ax1)) { 
      return false;
    }
    ax2.start = with.p1.AsPoint(1);
    ax2.direction = with.p1.Sub(with.p1, with.p2).Unity2D().AsVector();
    if (OnSide2D(ax2)) { 
      return false;
    }
    if (ax1.direction.Parallel2D(ax2.direction))
    {
      Range r1 = Project2D(ax1.direction);
      Range r2 = with.Project2D(ax2.direction);
      return r1.Overrlapping(r2);
    } else
    {
      return true;
    }
  }
  boolean Collide2D(Segment with, double halfThickness)
  {
    Line ax1 = new Line();
    Line ax2 = new Line();
    ax1.start = p1.AsPoint(1);
    ax1.direction = p1.Sub(p1, p2).Unity2D().AsVector();
    if (with.OnSide2D(ax1)) { 
      return false;
    }
    ax2.start = with.p1.AsPoint(1);
    ax2.direction = with.p1.Sub(with.p1, with.p2).Unity2D().AsVector();
    if (OnSide2D(ax2)) { 
      return false;
    }
    if (ax1.direction.Parallel2D(ax2.direction))
    {
      Range r1 = Project2D(ax1.direction);
      Range r2 = with.Project2D(ax2.direction);
      return r1.Overrlapping(r2, halfThickness);
    } else
    {
      return true;
    }
  }
  boolean Collide3D(Segment with)
  {
    Line ax1 = new Line();
    Line ax2 = new Line();
    ax1.start = p1.AsPoint(1);
    ax1.direction = p1.Sub(p1, p2).Unity2D().AsVector();
    if (with.OnSide2D(ax1)) { 
      return false;
    }
    ax2.start = with.p1.AsPoint(1);
    ax2.direction = with.p1.Sub(with.p1, with.p2).Unity2D().AsVector();
    if (OnSide2D(ax2)) { 
      return false;
    }
    if (ax1.direction.Parallel2D(ax2.direction))
    {
      Range r1 = Project2D(ax1.direction);
      Range r2 = with.Project2D(ax2.direction);
      boolean overlap = r1.Overrlapping(r2);
      if (overlap)
      {
        Segment szx = new Segment(new Vector(p1.x, p1.z), new Vector(p2.x, p2.z));
        Segment wszx = new Segment(new Vector(with.p1.x, with.p1.z), new Vector(with.p2.x, with.p2.z));
        return szx.Collide2D(wszx);
      }
      return overlap;
    } else
    {
      Segment szx = new Segment(new Vector(p1.x, p1.z), new Vector(p2.x, p2.z));
      Segment wszx = new Segment(new Vector(with.p1.x, with.p1.z), new Vector(with.p2.x, with.p2.z));
      return szx.Collide2D(wszx);
    }
  }
  boolean Collide3D(Segment with, double halfThickness)
  {
    Line ax1 = new Line();
    Line ax2 = new Line();
    ax1.start = p1.AsPoint(1);
    ax1.direction = p1.Sub(p1, p2).Unity2D().AsVector();
    if (with.OnSide2D(ax1)) { 
      return false;
    }
    ax2.start = with.p1.AsPoint(1);
    ax2.direction = with.p1.Sub(with.p1, with.p2).Unity2D().AsVector();
    if (OnSide2D(ax2)) { 
      return false;
    }
    if (ax1.direction.Parallel2D(ax2.direction))
    {
      Range r1 = Project2D(ax1.direction);
      Range r2 = with.Project2D(ax2.direction);
      boolean overlap = r1.Overrlapping(r2, halfThickness);
      if (overlap)
      {
        Segment szx = new Segment(new Vector(p1.x, p1.z), new Vector(p2.x, p2.z));
        Segment wszx = new Segment(new Vector(with.p1.x, with.p1.z), new Vector(with.p2.x, with.p2.z));
        return szx.Collide2D(wszx);
      }
      return overlap;
    } else
    {
      Segment szx = new Segment(new Vector(p1.x, p1.z), new Vector(p2.x, p2.z));
      Segment wszx = new Segment(new Vector(with.p1.x, with.p1.z), new Vector(with.p2.x, with.p2.z));
      return szx.Collide2D(wszx, halfThickness);
    }
  }
  boolean Collide2D(Line line)
  {
    return !OnSide2D(line);
  }
  boolean Collide3D(Line line)
  {
    return !OnSide3D(line);
  }
  double Slope2D(byte axis)
  {
    Cartesian p = p1.Sub(p1, p2);
    if (p.x != 0)
    {
      if (axis == AxisX)
        return (p.y / p.x) * p.y;
      else if (axis == AxisY)
        return (p.y) + p.Length2D();
    } else
      if (axis == AxisY)
        return (p.y) + p.Length2D();
    return 1;
  }
  Vector Parametric2D(double t)
  {
    return p1.Add(p1, p2.Unity2D().Scaled( t)).AsVector();
  }
  Vector Parametric3D(double t)
  {
    return p1.Add(p1, p2.Unity().Scaled(t)).AsVector();
  }
  Point ParametricIntersect2D(double t1, Segment  s2, double t2)
  {
    Cartesian v1 = p1.Add(p1, p2.Unity2D().Scaled(t1));
    Cartesian v2 = s2.p1.Add(s2.p1, s2.p2.Unity2D().Scaled( t2));
    Point solve = new Point(v1.x * t1 - v2.x * t2, v1.y * t1 + v2.y * t2);
    return solve;
  }
  Point ParametricIntersect3D(double t1, Segment s2, double t2)
  {
    Cartesian v1 = p1.Add(p1, p2.Unity().Scaled(t1));
    Cartesian v2 = s2.p1.Add(s2.p1, s2.p2.Unity().Scaled(t2));
    Point solve = new Point(v1.x * t1 - v2.x * t2, v1.y * t1 + v2.y * t2, v1.z * t1 - v2.z * t2);
    return solve;
  }
  Segment Add(Segment seg, Vector offset)
  {
    return new Segment(seg.p1.Add(seg.p1, offset).AsVector(), seg.p2.Add(seg.p2, offset).AsVector());
  }
  Segment Sub(Segment seg, Vector offset)
  {
    return new Segment(seg.p1.Sub(seg.p1, offset).AsVector(), seg.p2.Sub(seg.p2, offset).AsVector());
  }
};


class Circle extends Geometrical
{
  Point center;
  double radius;
  Circle(Point at, double rad)
  {
    center = at;
    radius = rad;
  }
  Circle()
  {

    center = new Point(0, 0);
    radius = 1;
  }
  boolean Collide2D(Circle with)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.Length2DSquared();
    return len < r2;
  }
  boolean Collide3D(Circle with)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.LengthSquared();
    return len < r2;
  }
  boolean Collide2D(Point p)
  {
    double r2 = radius * radius;
    Vector d = center.Sub(center, p);
    double dis2 = d.Length2DSquared();
    return dis2 < r2;
  }
  boolean Collide3D(Point p)
  {
    double r2 = radius * radius;
    Vector d = center.Sub(center, p);
    double dis2 = d.LengthSquared();
    return dis2 < r2;
  }
  boolean Collide2D(Point p, double r)
  {
    double r2 = radius  +  r;
    Vector d = center.Sub(center, p);
    double dis2 = d.Length2D();
    return dis2 < r2;
  }
  boolean Collide3D(Point p, double r)
  {
    double r2 = radius * radius + r * r;
    Vector d = center.Sub(center, p);
    double dis2 = d.LengthSquared();
    return dis2 < r2;
  }

  boolean Collide2D(Line ln)
  {
    Vector d =center.Sub( center, ln.start);
    Vector p = d.Project2D(ln.direction);
    Point nearest = ln.start.Add(ln.start, p).AsPoint(ln.start.r);
    return Collide2D(nearest);
  }
  boolean Collide3D(Line ln)
  {
    Vector d = center.Sub(center, ln.start);
    Vector p = d.Project(ln.direction);
    Point nearest = ln.start.Add(ln.start, p).AsPoint(ln.start.r);
    return Collide3D(nearest);
  }
  boolean Collide2D(Segment seg)
  {
    if (Collide2D(seg.p1.AsPoint(1))) return true;
    if (Collide2D(seg.p2.AsPoint(1))) return true;
    Cartesian d = seg.p1.Sub(seg.p1, seg.p2);
    Cartesian ln =center.Sub( center, seg.p1.AsPoint(1));
    Cartesian p = ln.Project2D(d);
    Point nearest = (seg.p1.Add(seg.p1, p)).AsPoint(1);
    return Collide2D(nearest) && p.Length2DSquared() <= d.Length2DSquared() && p.DotProduct2D(p, d) >= 0;
  }
  boolean Collide3D(Segment seg)
  {
    if (Collide3D(seg.p1.AsPoint(1))) return true;
    if (Collide3D(seg.p2.AsPoint(1))) return true;
    Cartesian d = seg.p1.Sub(seg.p1, seg.p2);
    Cartesian ln =center.Sub( center, seg.p1.AsPoint(1));
    Cartesian p = ln.Project(d);
    Point nearest = (seg.p1.Add(seg.p1, p)).AsPoint(1);
    return Collide3D(nearest) && p.LengthSquared() <= d.LengthSquared() && p.DotProduct(p, d) >= 0;
  }

  Circle Add(Circle c, Vector offset)
  {
    return new Circle(c.center.Add(c.center, offset).AsPoint(c.center.r), c.radius);
  }

  Circle Add(Circle c, double offset)
  {
    return new Circle(c.center, Math.abs(c.radius + offset));
  }
  Circle Sub(Circle c, double offset)
  {
    return new Circle(c.center, Math.abs(c.radius - offset));
  }
  Circle Sub(Circle c, Vector offset)
  {
    return new Circle(c.center.Sub(c.center, offset).AsPoint(c.center.r), c.radius);
  }
};

class Sphere extends Circle
{
  Normal pole;
  SphericalCoord direction;
  Sphere(Point at, double rad)
  {
    super(at, rad);

    pole = new Normal(0, 1, 0);
    direction = new SphericalCoord(new Point(0, 0, 1), AxisY);
  }
  Sphere()
  {
    pole = new Normal(0, 1, 0);
    direction = new SphericalCoord(new Point(0, 0, 1), AxisY);
  }
};

class Disk extends Circle
{
  double thickness;
  Normal planeNormal;
  //
  Disk(Point at, double rad, Normal plane)
  {
    super(at, rad);
    planeNormal = plane;
    thickness = Epsilon;
  }
  Disk(Point at, double rad, Normal plane, double thicken)

  {
    super(at, rad);
    planeNormal = plane;
    thickness = thicken;
  }
  Disk()
  {
    planeNormal = new Normal(0, 1, 0);
    thickness = Epsilon;
  }

  boolean Collide3D(Circle with)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.LengthSquared();
    if (len < r2)
    {
      double added = thickness + Epsilon;
      double d1 = dif.DotProduct(dif, planeNormal);
      double d2 = dif.DotProduct(dif, new Normal(0, 1, 0));
      double d3 = Math.abs(d1 + d2);
      return (d3 < added);
    } else
    {
      return false;
    }
  }
  boolean CollideWithSphere3D(Circle with)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.LengthSquared();
    if (len < r2)
    {
      double added = thickness + Epsilon;
      double d1 = dif.DotProduct(dif, planeNormal);
      return (d1 > -added && d1 < added);
    } else
    {
      return false;
    }
  }
  boolean Collide3D(Circle with, double thickness)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.LengthSquared();
    if (len < r2)
    {
      double added = this.thickness + thickness;
      double d1 = dif.DotProduct(dif, planeNormal);
      double d2 = dif.DotProduct(dif, new Normal(0, 1, 0));
      double d3 = Math.abs(d1 + d2);
      return (d3 < added);
    } else
    {
      return false;
    }
  }
  boolean Collide3D(Circle with, Normal norm, double halfThickness)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.LengthSquared();
    if (len < r2)
    {
      double added = this.thickness + halfThickness;
      double d1 = dif.DotProduct(dif, planeNormal);
      double d2 = dif.DotProduct(dif, norm);
      double d3 = Math.abs(d1 + d2);
      return (d3 < added);
    } else
    {
      return false;
    }
  }
  boolean Collide3D(Disk with)
  {
    double r2 = radius * radius + with.radius * with.radius;
    Vector dif = center.Sub(center, with.center);
    double len = dif.LengthSquared();
    if (len < r2)
    {
      double added = thickness + with.thickness;
      double d1 = dif.DotProduct(dif, planeNormal);
      double d2 = dif.DotProduct(dif, with.planeNormal);
      double d3 = Math.abs(d1 + d2);
      return (d3 < added);
    } else
    {
      return false;
    }
  }
};

class Box extends Geometrical
{
  Point origin;
  Vector size;

  Box(Point from, Vector sz)
  {
    origin = from;
    size = sz.Absolute();
  }
  Box()
  {
    origin = new Point(-0.5f, -0.5f);
    size = new Vector(1, 1);
  }
  Box Enlarge2D(Vector with)
  {
    Box e = new Box();
    e.origin.x = Math.min(origin.x, with.x);
    e.origin.y = Math.min(origin.y, with.y);
    e.size.x = Math.max(origin.x + size.x, with.x);
    e.size.y = Math.max(origin.y + size.y, with.y);
    e.size =  e.size.Sub(e.size, e.origin).AsVector();
    return e;
  }
  Box Enlarge3D(Vector with)
  {
    Box e = new Box();
    e.origin.x = Math.min(origin.x, with.x);
    e.origin.y = Math.min(origin.y, with.y);
    e.origin.z = Math.min(origin.z, with.z);
    e.size.x = Math.max(origin.x + size.x, with.x);
    e.size.y = Math.max(origin.y + size.y, with.y);
    e.size.z = Math.max(origin.z + size.z, with.z);
    e.size =e.size.Sub(e.size, e.origin).AsVector();
    return e;
  }
  Vector Corner2D(byte face, byte side)
  {
    Vector corner = origin.AsVector();
    switch (face)
    {

    case SideBack:
      {
        corner.y += size.y;
        break;
      }
    case SideFront:
      {
        break;
      }
    }
    switch (side)
    {
    case SideLeft:
      {
        break;
      }
    case SideRight:
      {
        corner.x += size.x;
        break;
      }
    }
    return corner;
  }
  Vector Corner3D(byte face, byte side, byte base)
  {
    Vector corner = Corner2D(face, side);
    switch (base)
    {
    case SideTop:
      {
        corner.z += size.z;
        break;
      }
    }
    return corner;
  }
  boolean SeparatingAxis2D(Segment axis)
  {
    Segment e1 = new Segment(); 
    Segment e2 = new Segment();
    Range r = new Range();
    Range r1 = new Range(); 
    Range r2 = new Range(); 
    Range pro = new Range();
    Vector n = axis.p1.Sub(axis.p1, axis.p2).AsVector();
    e1.p1 = Corner2D(SideFront, SideLeft);
    e1.p2 = Corner2D(SideBack, SideRight);
    e2.p1 = Corner2D(SideBack, SideLeft);
    e2.p2 = Corner2D(SideFront, SideRight);
    r1 = e1.Project2D(n);
    r2 = e2.Project2D(n);
    pro = new Range(r1, r2);
    r = axis.Project2D(n);
    return !r.Overrlapping(pro);
  }
  boolean SeparatingAxis3D(Segment axis)
  {
    Segment e1 = new Segment(); 
    Segment e2 = new Segment(); 
    Segment e3 = new Segment();
    Range r = new Range();
    Range r1 = new Range(); 
    Range r2 = new Range(); 
    Range r3 = new Range(); 
    Range pro = new Range();
    Vector n = axis.p1.Sub(axis.p1, axis.p2).AsVector();
    e1.p1 = Corner3D(SideFront, SideLeft, SideTop);
    e1.p2 = Corner3D(SideBack, SideRight, SideBottom);
    e2.p1 = Corner3D(SideBack, SideLeft, SideBottom);
    e2.p2 = Corner3D(SideFront, SideRight, SideTop);
    e3.p1 = Corner3D(SideFront, SideLeft, SideBottom);
    e3.p2 = Corner3D(SideBack, SideRight, SideTop);
    r1 = e1.Project3D(n);
    r2 = e2.Project3D(n);
    r3 = e3.Project3D(n);
    pro = new Range(new Range(r1, r2), r3);
    r = axis.Project2D(n);
    return !r.Overrlapping(pro);
  }
  boolean Collide2D(Box with)
  {
    Point extx = new Point(origin.x, origin.x + size.x);
    Point wextx = new Point(with.origin.x, with.origin.x + with.size.x);
    Point exty = new Point(origin.y, origin.y + size.y);
    Point wexty = new Point(with.origin.y, with.origin.y + with.size.y);
    return Overlapping2D(extx, wextx) && Overlapping2D(exty, wexty);
  }
  boolean Collide3D(Box with)
  {
    Point extx = new Point(origin.x, origin.x + size.x);
    Point wextx = new Point(with.origin.x, with.origin.x + with.size.x);
    Point exty = new Point(origin.y, origin.y + size.y);
    Point wexty = new Point(with.origin.y, with.origin.y + with.size.y);
    Point extz = new Point(origin.z, origin.z + size.z);
    Point wextz = new Point(with.origin.z, with.origin.z + with.size.z);
    return Overlapping2D(extx, wextx) && Overlapping2D(exty, wexty) && Overlapping2D(extz, wextz);
  }
  Point Clamp2D(Point p)
  {
    Range x = new Range(origin.x, origin.x + size.x);
    Range y = new Range(origin.y, origin.y + size.y);
    return new Point(x.Clamp(p.x), y.Clamp(p.y));
  }
  Point Clamp3D(Point p)
  {
    Range x = new Range(origin.x, origin.x + size.x);
    Range y = new Range(origin.y, origin.y + size.y);
    Range z = new Range(origin.z, origin.z + size.z);
    return new Point(x.Clamp(p.x), y.Clamp(p.y), z.Clamp(p.z));
  }
  boolean Collide2D(Point p)
  {
    double left = origin.x;
    double right = origin.x + size.x;
    double back = origin.y;
    double front = origin.y + size.y;

    if (p.x < left) { 
      return false;
    } else
      if (p.x > right) { 
        return false;
      } else
        if (p.y > front) { 
          return false;
        } else
          if (p.y <back) { 
            return false;
          }
    //System.Console.WriteLine("TEST" + p.x + " " + p.y + " " + left + " " + right + " " + back + " " + front);
    return true;
  }
  boolean Collide3D(Point p)
  {
    double left = origin.x;
    double right = origin.x + size.x;
    double back = origin.y;
    double front = origin.y + size.y;
    double bottom = origin.z;
    double top = origin.z + size.z;

    if (p.x < left) { 
      return false;
    } else
      if (p.x > right) { 
        return false;
      } else
        if (p.y > front) { 
          return false;
        } else
          if (p.y < back) { 
            return false;
          } else
            if (p.z > top) { 
              return false;
            } else
              if (p.z < bottom) { 
                return false;
              } else
                return true;
  }
  boolean Collide2D(Line ln)
  {
    Cartesian n = ln.direction.Rotate90Degrees2D();
    double dp1, dp2, dp3, dp4;
    Point c1 = origin;
    Point c2 = origin.Add(origin, size).AsPoint(origin.r);
    Point c3 = new Point(c2.x, c1.y);
    Point c4 = new Point(c1.x, c2.y);
    Vector v1 = c1.Sub(c1, ln.start);
    Vector v2 = c2.Sub(c2, ln.start);
    Vector v3 = c3.Sub(c3, ln.start);
    Vector v4 = c4.Sub(c4, ln.start);
    dp1 = n.DotProduct2D(n, v1);
    dp2 = n.DotProduct2D(n, v2);
    dp3 = n.DotProduct2D(n, v3);
    dp4 = n.DotProduct2D(n, v4);
    return (dp1 * dp2 <= 0) || (dp2 * dp3 <= 0) || (dp3 * dp4 <= 0);
  }
  boolean Collide3D(Line ln)
  {
    Vector n = ln.start.CrossProduct(ln.start, ln.direction).Unity().AsVector();
    double dp1, dp2, dp3, dp4, dp5, dp6;
    Point c1 = origin;
    Point c2 = origin.Add(origin, size).AsPoint(origin.r);
    Point c3 = new Point(c2.x, c1.y, c1.z);
    Point c4 = new Point(c1.x, c2.y, c1.z);
    Point c5 = new Point(c2.x, c1.y, c2.z);
    Point c6 = new Point(c1.x, c2.y, c2.z);
    Vector v1 = c1.Sub(c1, ln.start);
    Vector v2 = c2.Sub(c2, ln.start);
    Vector v3 = c3.Sub(c3, ln.start);
    Vector v4 = c4.Sub(c4, ln.start);
    Vector v5 = c5.Sub(c5, ln.start);
    Vector v6 = c6.Sub(c6, ln.start);
    dp1 = n.DotProduct(n, v1);
    dp2 = n.DotProduct(n, v2);
    dp3 = n.DotProduct(n, v3);
    dp4 = n.DotProduct(n, v4);
    dp5 = n.DotProduct(n, v5);
    dp6 = n.DotProduct(n, v6);
    return (dp1 * dp2 <= 0) || (dp2 * dp3 <= 0) || (dp3 * dp4 <= 0) || (dp4 * dp5 <= 0) || (dp5 * dp6 <= 0);
  }
  boolean Collide2D(Segment seg)
  {
    Line line = new Line();
    line.start = seg.p1.AsPoint(1);
    line.direction = seg.p2.Sub(seg.p2, seg.p1).Unity2D().AsVector();
    if (!Collide2D(line))
    {
      return false;
    }
    Range r = new Range();
    Range s = new Range();
    r.minim = origin.x;
    r.maxim = origin.x + size.x;
    s.minim = seg.p1.x;
    s.maxim = seg.p2.x;
    s = new Range(r, s);
    if (!r.Overrlapping(s)) { 
      return false;
    }
    r.minim = origin.y;
    r.maxim = origin.y + size.y;
    s.minim = seg.p1.y;
    s.maxim = seg.p2.y;
    s = new Range(r, s);
    if (!r.Overrlapping(s)) { 
      return false;
    }
    return true;
  }
  boolean Collide3D(Segment seg)
  {
    Line line = new Line();
    line.start = seg.p1.AsPoint(1);
    line.direction = seg.p2.Sub(seg.p2, seg.p1).Unity().AsVector();
    if (!Collide3D(line))
    {
      return false;
    }
    Range r = new Range();
    Range s = new Range();
    r.minim = origin.x;
    r.maxim = origin.x + size.x;
    s.minim = seg.p1.x;
    s.maxim = seg.p2.x;
    s = new Range(r, s);
    if (!r.Overrlapping(s)) { 
      return false;
    }
    r.minim = origin.y;
    r.maxim = origin.y + size.y;
    s.minim = seg.p1.y;
    s.maxim = seg.p2.y;
    s = new Range(r, s);
    if (!r.Overrlapping(s)) { 
      return false;
    }
    r.minim = origin.z;
    r.maxim = origin.z + size.z;
    s.minim = seg.p1.z;
    s.maxim = seg.p2.z;
    s = new Range(r, s);
    if (!r.Overrlapping(s)) { 
      return false;
    }
    return true;
  }
  boolean Collide2D(Circle circ)
  {
    Point clamp = Clamp2D(circ.center);
    return circ.Collide2D(clamp);
  }
  boolean Collide3D(Circle circ)
  {
    Point clamp = Clamp3D(circ.center);
    return circ.Collide3D(clamp);
  }
  Box Add(Box p, Vector offset)
  {
    Point at=p.origin.Add(p.origin, offset).AsPoint(p.origin.r);
    return new Box(at, p.size);
  }
  Box Sub(Box p, Point offset)
  {
    Cartesian at = p.origin.Sub(p.origin, offset);

    return new Box(at.AsPoint(p.origin.r), p.size);
  }
};

class Rectalinear extends Circle
{
  Vector halfExtent;
  SphericalCoord angle;
  Rectalinear(Point at, Vector hext)
  {
    super(at, hext.Length());
    halfExtent = hext;
    angle = new SphericalCoord(at.Add(at, hext), AxisZ);
  }
  Rectalinear()
  {
    halfExtent = new Vector(0.5f, 0.5f);
    angle = new SphericalCoord(halfExtent, AxisZ);
  }
  Vector HalfCorner2D(byte face, byte side)
  {

    Vector corner = halfExtent;
    switch (face)
    {

    case SideBack:
      {
        corner.y -= corner.y;
        break;
      }
    case SideFront:
      {
        corner.y = -corner.y;
        break;
      }
    }
    switch (side)
    {
    case SideLeft:
      {

        corner.x = -corner.x;
        break;
      }
    case SideRight:
      {
        corner.x -= corner.x;
        break;
      }
    }
    return corner;
  }
  Vector HalfCorner3D(byte face, byte side, byte base)
  {

    Vector corner = HalfCorner2D(face, side);
    switch (base)
    {

    case SideTop:
      {
        corner.z -= corner.z;
        break;
      }
    case SideBottom:
      {
        corner.z = -corner.z;
        break;
      }
    }
    return corner;
  }
  Box Hull2D()
  {
    Box h = new Box(center, new Vector(0, 0));
    Vector c;
    c = HalfCorner2D(SideFront, SideLeft).Rotate2D(angle.theta).AsVector();
    h.Enlarge2D(c);
    c = HalfCorner2D(SideBack, SideRight).Rotate2D(angle.theta).AsVector();
    h.Enlarge2D(c);
    c = HalfCorner2D(SideFront, SideRight).Rotate2D(angle.theta).AsVector();
    h.Enlarge2D(c);
    c = HalfCorner2D(SideBack, SideLeft).Rotate2D(angle.theta).AsVector();
    h.Enlarge2D(c);
    return h;
  }
  Box Hull3D()
  {
    Box h = new Box(center, new Vector(0, 0));
    Vector c;
    c = HalfCorner3D(SideFront, SideLeft, SideTop);
    SphericalCoord ang = new SphericalCoord(c, AxisZ);
    ang.phi += angle.phi;
    ang.theta += angle.theta;
    h.Enlarge3D(ang.AtZPole());

    c = HalfCorner3D(SideBack, SideRight, SideBottom);
    ang = new SphericalCoord(c, AxisZ);
    ang.phi += angle.phi;
    ang.theta += angle.theta;

    h.Enlarge3D(ang.AtZPole());
    c = HalfCorner3D(SideFront, SideRight, SideTop);
    ang = new SphericalCoord(c, AxisZ);
    ang.phi += angle.phi;
    ang.theta += angle.theta;

    h.Enlarge3D(ang.AtZPole());
    c = HalfCorner3D(SideBack, SideLeft, SideBottom);
    ang = new SphericalCoord(c, AxisZ);
    ang.phi += angle.phi;
    ang.theta += angle.theta;

    h.Enlarge3D(ang.AtZPole());
    c = HalfCorner3D(SideBack, SideRight, SideTop);
    ang = new SphericalCoord(c, AxisZ);
    ang.phi += angle.phi;
    ang.theta += angle.theta;

    h.Enlarge3D(ang.AtZPole());
    c = HalfCorner3D(SideFront, SideLeft, SideBottom);
    ang = new SphericalCoord(c, AxisZ);
    ang.phi += angle.phi;
    ang.theta += angle.theta;

    h.Enlarge3D(ang.AtZPole());
    return h;
  }
  Range Hull2D(Range r1, Range r2)
  {
    return new Range(r1, r2);
  }
  Range Hull3D(Range r1, Range r2, Range r3)
  {
    Range p1 = new Range(r1, r2);
    Range p2 = new Range(p1, r3);
    return p2;
  }
  boolean Collide2D(Line line)
  {
    if (super.Collide2D(line))
    {
      Box r = new Box();
      r.origin.x = 0;
      r.origin.y = 0;
      r.size = halfExtent.Scaled(2);
      Line l2 = new Line();
      double rad=l2.start.r;
      l2.start = line.start.Sub(line.start, center).AsPoint(rad);
      l2.start = l2.start.Rotate2D(-angle.theta).AsPoint(rad);
      l2.start = l2.start.Add(l2.start, halfExtent).AsPoint(rad);
      l2.direction = line.direction.Rotate2D(-angle.theta).AsVector();
      return r.Collide2D(l2);
    }
    return false;
  }
  boolean Collide2D(Segment seg)
  {
    if (super.Collide2D(seg))
    {

      Box r = new Box();
      r.origin.x = 0;
      r.origin.y = 0;
      r.size = halfExtent.Scaled(2);
      Segment s2 = new Segment();
      s2.p1 = seg.p1.Sub(seg.p1, center).AsVector();
      s2.p1 = s2.p1.Rotate2D(-angle.theta).AsVector();
      s2.p1 = s2.p1.Add(s2.p1, halfExtent).AsVector();
      s2.p2 = seg.p2.Sub(seg.p2, center).AsVector();
      s2.p2 = s2.p2.Rotate2D(-angle.theta).AsVector();
      s2.p2 = s2.p2.Add(s2.p2, halfExtent).AsVector();
      return r.Collide2D(s2);
    }
    return false;
  }
  boolean Collide3D(Segment seg)
  {
    if (super.Collide3D(seg))
    {

      Box r = new Box();
      r.origin.x = 0;
      r.origin.y = 0;
      r.size = halfExtent.Scaled(2);
      Segment s2 = new Segment();
      s2.p1 = seg.p1.Sub(seg.p1, center).AsVector();
      s2.p1 = s2.p1.Rotate3D(-angle.phi, -angle.theta).AsVector();
      s2.p1 = s2.p1.Add(s2.p1, halfExtent).AsVector();
      s2.p2 = seg.p2.Sub(seg.p2, center).AsVector();
      s2.p2 = s2.p2.Rotate3D(-angle.phi, -angle.theta).AsVector();
      s2.p2 = s2.p2.Add(s2.p2, halfExtent).AsVector();
      return r.Collide3D(s2);
    }
    return false;
  }
  boolean Collide3D(Line line)
  {
    if (super.Collide3D(line))
    {

      Box r = new Box();
      r.origin.x = 0;
      r.origin.y = 0;
      r.size = halfExtent.Scaled(2);
      Line l2 = new Line();
      double rad=l2.start.r;
      l2.start = line.start.Add(line.start, center).AsPoint(rad);
      l2.start = l2.start.Rotate3D(-angle.phi, -angle.theta).AsPoint(rad);
      l2.start = l2.start.Add(l2.start, halfExtent).AsPoint(rad);
      l2.direction = line.direction.Rotate3D(-angle.phi, -angle.theta).AsVector();
      return r.Collide3D(l2);
    }
    return false;
  }
  boolean Collide2D(Point with)
  {
    if (super.Collide2D(with))
    {

      Box r = new Box();
      r.origin.x = 0;
      r.origin.y = 0;
      r.size = halfExtent.Scaled(2);
      Vector p = with.Sub(with, center);
      p = p.Rotate2D(-angle.theta).AsVector();
      p = p.Add(p, halfExtent).AsVector();
      return r.Collide2D(p.AsPoint(with.r));
    }
    return false;
  }
  boolean Collide3D(Point with)
  {
    if (super.Collide3D(with))
    {

      Box r = new Box();
      r.origin.x = 0;
      r.origin.y = 0;
      r.size = halfExtent.Scaled(2);
      Vector p = with.Sub(with, center).AsVector();
      p = p.Rotate3D(angle.phi, angle.theta).AsVector();
      p = p.Add(p, halfExtent).AsVector();
      return r.Collide3D(p.AsPoint(with.r));
    }
    return false;
  }
  Segment Edge2D(byte side)
  {
    Segment edge = new Segment();
    Vector e1 = new Vector(halfExtent.x, halfExtent.y);
    Vector e2 = new Vector(halfExtent.x, halfExtent.y);
    switch (side)
    {
    case SideFront: 
      { 
        e1.x = -e1.x; 
        break;
      }
    case SideRight: 
      { 
        e2.y = -e2.y; 
        break;
      }
    case SideLeft: 
      { 
        e1.y = -e1.y; 
        e2 = e2.Negative(); 
        break;
      }
    case SideBack: 
      { 
        e1 = e1.Negative(); 
        e2.x = -e2.x; 
        break;
      }
    }
    e1 = e1.Rotate2D(angle.theta).AsVector();
    e1 = e1.Add(e1, center).AsVector();
    e2 = e2.Rotate2D(angle.theta).AsVector();
    e2 =e2.Add(e2, center).AsVector();
    edge.p1 = e1;
    edge.p2 = e2;
    return edge;
  }
  Segment Edge3D(Point at, byte side)
  {
    Segment edge = new Segment();
    Vector e1 = new Vector(at.x, at.y);
    Vector e2 = new Vector(at.x, at.y);
    switch (side)
    {
    case SideFront: 
      { 
        e1.x = -e1.x; 
        break;
      }
    case SideRight: 
      { 
        e2.y = -e2.y; 
        break;
      }
    case SideTop: 
      { 
        e2.z = -e2.z; 
        break;
      }
    case SideBottom: 
      { 
        e1 = e1.Negative(); 
        e2.z = -e2.z; 
        break;
      }
    case SideLeft: 
      { 
        e1.y = -e1.y; 
        e2 = e2.Negative(); 
        break;
      }
    case SideBack: 
      { 
        e1 = e1.Negative(); 
        e2.x = -e2.x; 
        break;
      }
    }
    e1 = e1.Add(e1, center).AsVector();
    e2 = e2.Add(e2, center).AsVector();
    edge.p1 = e1;
    edge.p2 = e2;
    return edge;
  }
  boolean SeparatingAxis2D(Segment axis)
  {
    Range ar = new Range(); 
    Range r0 = new Range(); 
    Range r2 = new Range(); 
    Range project = new Range();
    Segment e0 = Edge2D(SideFront);
    Segment e2 = Edge2D(SideLeft);
    Vector n = axis.p1.Sub(axis.p1, axis.p2).AsVector();
    ar = axis.Project2D(n);
    r0 = e0.Project2D(n);
    r2 = e2.Project2D(n);
    project = Hull2D(r0, r2);
    return !ar.Overrlapping(project);
  }
  boolean SeparatingAxis3D(Segment axis)
  {
    Range ar = new Range(); 
    Range r0 = new Range(); 
    Range r1 = new Range(); 
    Range r2 = new Range(); 
    Range pro = new Range();

    Segment e0 = Edge2D(SideFront);
    Segment e1 = Edge2D(SideLeft);
    Segment e2 = Edge2D(SideTop);
    Vector n = axis.p1.Sub(axis.p1, axis.p2).AsVector();
    ar = axis.Project3D(n);
    r0 = e0.Project3D(n);
    r1 = e1.Project3D(n);
    r2 = e2.Project3D(n);
    pro = Hull3D(r0, r1, r2);
    return !ar.Overrlapping(pro);
  }
  boolean Collide2D(Rectalinear with)
  {
    Segment edge = Edge2D(SideFront);
    if (with.SeparatingAxis2D(edge)) { 
      return false;
    }
    edge = Edge2D(SideLeft);
    if (with.SeparatingAxis2D(edge)) { 
      return false;
    }
    edge = Edge2D(SideRight);
    if (with.SeparatingAxis2D(edge)) { 
      return false;
    }
    edge = Edge2D(SideBack);
    return !with.SeparatingAxis2D(edge);
  }
  boolean Collide3D(Rectalinear with)
  {
    Range cz = new Range(center.z - halfExtent.z, center.z + halfExtent.z);
    Range czw = new Range(with.center.z - with.halfExtent.z, with.center.z + with.halfExtent.z);
    if (cz.Overrlapping(czw))
    {

      angle.p = radius;
      Point at = angle.AtZPole().AsPoint(1);
      Segment edge = Edge3D(at, SideFront);
      if (with.SeparatingAxis3D(edge)) { 
        return false;
      }
      edge = Edge3D(at, SideLeft);
      if (with.SeparatingAxis3D(edge)) { 
        return false;
      }
      edge = Edge3D(at, SideTop);
      if (with.SeparatingAxis3D(edge)) { 
        return false;
      }
      edge = Edge3D(at, SideRight);
      if (with.SeparatingAxis3D(edge)) { 
        return false;
      }
      edge = Edge3D(at, SideBack);
      if (with.SeparatingAxis3D(edge)) { 
        return false;
      }
      edge = Edge3D(at, SideBottom);
      if (with.SeparatingAxis3D(edge)) { 
        return false;
      }

      return true;
    }

    return false;
  }

  boolean CollideWith2D(Circle with)
  {
    if (super.Collide2D(with))
    {
      Box r = new Box();
      r.size = new Vector(halfExtent.x * 2, halfExtent.y * 2);
      Circle c = new Circle();
      c.radius = with.radius;
      Vector d = with.center.Sub(with.center, center);
      d = d.Rotate2D(-angle.theta).AsVector();
      c.center = (d.Add(d, halfExtent)).AsPoint(c.radius);
      return r.Collide2D(c);
    }
    return false;
  }
  boolean CollideWith3D(Circle with)
  {
    if (super.Collide3D(with))
    {
      Box r = new Box();
      r.size = new Vector(halfExtent.x * 2, halfExtent.y * 2, halfExtent.z * 2);
      Circle c = new Circle();
      c.radius = with.radius;
      Vector d = with.center.Sub(with.center, center);
      Vector at = angle.AtZPole();
      c.center = d.Sub(d, at).AsPoint(c.center.r);
      c.center=c.center.Add(c.center, halfExtent).AsPoint(c.radius);
      return r.Collide2D(c);
    }
    return false;
  }
  boolean CollideWith2D(Rectalinear with)
  {
    if (super.Collide2D(with))
    {
      return Collide2D(with);
    }
    return false;
  }
  boolean CollideWith3D(Rectalinear with)
  {
    if (super.Collide3D(with))
    {
      return Collide3D(with);
    }
    return false;
  }
};

class BoundingBox extends Geometrical
{
  Vector pMin;
  Vector pMax;
  BoundingBox()
  {
    pMin = new Vector(MAX_MAP, MAX_MAP, MAX_MAP);
    pMax = new Vector(-MAX_MAP, -MAX_MAP, -MAX_MAP);
  }
  BoundingBox(Point p)
  {
    pMin = new Vector(p.x, p.y, p.z, p.r);
    pMax = new Vector(p.x, p.y, p.z, p.r);
  }
  BoundingBox(Point p1, Point p2)
  {
    pMin = new Vector(Math.min(p1.x, p2.x), 
      Math.min(p1.y, p2.y), 
      Math.min(p1.z, p2.z));
    pMax = new Vector(Math.max(p1.x, p2.x), 
      Math.max(p1.y, p2.y), 
      Math.max(p1.z, p2.z));
  }
  BoundingBox(double mx, double my, double mz)
  {
    pMin = new Vector(-mx, -my, -mz);
    pMax = new Vector(mx, my, mz);
  }
  BoundingBox(double mix, double miy, double miz, double max, double may, double maz)
  {
    pMin = new Vector(-mix, -miy, -miz);
    pMax = new Vector(max, may, maz);
  }
  boolean Overlaps(BoundingBox b)
  {
    boolean x = (pMax.x >= b.pMin.x) && (pMin.x <= b.pMax.x);
    boolean y = (pMax.y >= b.pMin.y) && (pMin.y <= b.pMax.y);
    boolean z = (pMax.z >= b.pMin.z) && (pMin.z <= b.pMax.z);
    return (x && y && z);
  }
  boolean Inside(Point pt)
  {
    return (pt.x >= pMin.x && pt.x <= pMax.x &&
      pt.y >= pMin.y && pt.y <= pMax.y &&
      pt.z >= pMin.z && pt.z <= pMax.z);
  }
  boolean Inside(Point pt, double r)
  {
    return (pt.x - r >= pMin.x && pt.x + r <= pMax.x &&
      pt.y - r >= pMin.y && pt.y + r <= pMax.y &&
      pt.z - r >= pMin.z && pt.z + r <= pMax.z);
  }
  void Expand(double delta)
  {
    pMin = new Vector(pMin.x-delta, pMin.y-delta, pMin.z-delta);
    pMax = new Vector(pMax.x+delta, pMax.y+delta, pMax.z+delta);
  }
  double Volume()
  {
    Vector d = pMax.Sub(pMax, pMin).AsVector();
    return d.x * d.y * d.z;
  }
  int MaximumExtent()
  {
    Vector diag = pMax.Sub(pMax, pMin).AsVector();
    if (diag.x > diag.y && diag.x > diag.z)
      return 0;
    else if (diag.y > diag.z)
      return 1;
    else
      return 2;
  }
}
