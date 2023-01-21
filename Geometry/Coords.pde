class CylindricalCoord
{

  double r;
  double xyz;
  double theta;
  CylindricalCoord()
  {
    xyz = 0;
    r = 1;
    theta = 0;
  }
  CylindricalCoord(double radius, double angle, double axis)
  {
    xyz = axis;
    r = radius;
    theta = angle;
  }
  CylindricalCoord(Cartesian at, byte poleAxis)
  {
    if (poleAxis == AxisZ)
    {
      r = at.Length2D();
      theta = (double)(Math.atan2(at.y, at.x));
      xyz = at.z;
    } else

      if (poleAxis == AxisY)
      {
        r = (double)(Math.sqrt(at.x*at.x+at.z*at.z));
        theta = (double)(Math.atan2(at.z, at.x));
        xyz = at.y;
      } else
      {
        r = (double)(Math.sqrt(at.y * at.y + at.z * at.z));
        theta = (double)(Math.atan2(at.z, at.y));
        xyz = at.x;
      }
  }
  Point AtZPole()
  {
    return new Point(r * (double)Math.cos(theta), r * (double)Math.sin(theta), xyz);
  }
  Point AtYPole()
  {
    return new Point(r * (double)Math.cos(theta), xyz, r * (double)Math.sin(theta));
  }
  Point AtXPole()
  {
    return new Point(xyz, r * (double)Math.cos(theta), r * (double)Math.sin(theta));
  }
  Normal ToZPole()
  {
    Normal n = new Normal(r * (double)Math.cos(theta), r * (double)Math.sin(theta), xyz);
    n.Normalize();
    return n;
  }
  Normal ToYPole()
  {
    Normal n = new Normal(r * (double)Math.cos(theta), xyz, r * (double)Math.sin(theta));
    n.Normalize();
    return n;
  }
  Normal ToXPole()
  {
    Normal n = new Normal(xyz, r * (double)Math.cos(theta), r * (double)Math.sin(theta));
    n.Normalize();
    return n;
  }
  double Sine(double invr)
  {
    return xyz * invr;
  }
  double Cosine(double invr)
  {
    return xyz * invr;
  }
}
class Plane
{
  double distance;
  Normal normal;
  Plane()
  {
    distance = 0;
    normal = new Normal(0, 1, 0);
  }
  Plane(Point pos, Normal norm)
  {
    normal = norm;
    distance = -(normal.DotProduct(normal, pos));
  }
  Plane(Normal dir, double d, Normal norm)
  {
    normal = norm;
    distance = -(normal.DotProduct(norm, dir.Scaled(d)));
  }
  Plane(double d, Normal norm)
  {

    normal = norm;
    distance = d;
  }
  Plane(Normal norm, double d)
  {

    normal = norm;
    distance = d;
  }
  Plane(Point pos1, Point pos2, Point pos3)
  {
    Normal n = pos1.CrossProduct(pos1.Sub(pos2, pos1), pos1.Sub(pos3, pos1)).Normal();
    normal = new Normal(n.x, n.y, n.z);
    distance = -(normal.DotProduct(n, pos1));
  }
  Plane(Vector pos1, Vector pos2, Vector pos3)
  {
    Normal n = pos1.CrossProduct(pos1.Sub(pos2, pos1), pos1.Sub(pos3, pos1)).Normal();
    normal = new Normal(n.x, n.y, n.z);
    distance = -(normal.DotProduct(n, pos1));
  }
  void Flip()
  {
    normal = normal.Negative().AsNormal();
  }
  double Product(Point p)
  {

    return normal.DotProduct( normal, p) + distance;
  }

  boolean Front(Point p, double product)
  {
    double dp = product;
    if (dp > E)
    {
      return true;
    }
    return false;
  }
  boolean Back(Point p, double product)
  {
    double dp = product;
    if (dp < -Epsilon)
    {
      return true;
    }
    return false;
  }

  boolean On(Point p, double product)
  {
    double dp = product;
    if (dp < E && dp > -E)
    {
      return true;
    }
    return false;
  }
  boolean Front(Point p)
  {
    double dp = normal.DotProduct(normal, p) + distance;
    if (dp > E)
    {
      return true;
    }
    return false;
  }
  boolean Back(Point p)
  {
    double dp = normal.DotProduct( normal, p) + distance;
    if (dp < -E)
    {
      return true;
    }
    return false;
  }

  boolean On(Point p)
  {
    double dp = normal.DotProduct( normal, p) + distance;
    if (dp < E && dp > -E)
    {
      return true;
    }
    return false;
  }

  boolean Front(Point p, double product, double thickness)
  {
    double dp = product;
    if (dp > thickness)
    {
      return true;
    }
    return false;
  }
  boolean Back(Point p, double product, double thickness)
  {
    double dp = product;
    if (dp < -thickness)
    {
      return true;
    }
    return false;
  }

  boolean On(Point p, double product, double thickness)
  {
    double dp = product;
    if (dp < thickness && dp > -thickness)
    {
      return true;
    }
    return false;
  }
}
class SphericalCoord
{

  double p;
  double theta;
  double phi;
  SphericalCoord()
  {
    p = 1;
    theta = 0;
    phi = 0;
  }
  SphericalCoord(SphericalCoord from, SphericalCoord to)
  {
    p = (from.p+to.p)*0.5f;
    theta = (from.theta+to.theta)*0.5f;
    phi = (from.phi + to.phi) * 0.5f;
  }
  SphericalCoord(double xang, double yang, double rad)
  {
    phi = yang;
    theta = xang;
    p = rad;
  }
  SphericalCoord(Cartesian at, byte pole)
  {

    if (pole == AxisZ)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.y * at.y);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.y, at.x);
      phi = (double)Math.acos(at.z / p);
    } else if (pole == AxisY)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.x);
      phi = (double)Math.acos(at.y / p);
    } else
    {
      double r = (double)Math.sqrt(at.y * at.y + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.y);
      phi = (double)Math.acos(at.x / p);
    }
  }
  SphericalCoord(Point at, byte pole)
  {

    if (pole == AxisZ)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.y * at.y);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.y, at.x);
      phi = (double)Math.acos(at.z / p);
    } else if (pole == AxisY)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.x);
      phi = (double)Math.acos(at.y / p);
    } else
    {
      double r = (double)Math.sqrt(at.y * at.y + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.y);
      phi = (double)Math.acos(at.x / p);
    }
  }
  SphericalCoord(Vector at, byte pole)
  {

    if (pole == AxisZ)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.y * at.y);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.y, at.x);
      phi = (double)Math.acos(at.z / p);
    } else if (pole == AxisY)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.x);
      phi = (double)Math.acos(at.y / p);
    } else
    {
      double r = (double)Math.sqrt(at.y * at.y + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.y);
      phi = (double)Math.acos(at.x / p);
    }
  }
  SphericalCoord(Normal at, byte pole)
  {

    if (pole == AxisZ)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.y * at.y);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.y, at.x);
      phi = (double)Math.acos(at.z / p);
    } else if (pole == AxisY)
    {
      double r = (double)Math.sqrt(at.x * at.x + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.x);
      phi = (double)Math.acos(at.y / p);
    } else
    {
      double r = (double)Math.sqrt(at.y * at.y + at.z * at.z);
      p = (double)at.Length();
      theta = (double)Math.atan2(at.z, at.y);
      phi = (double)Math.acos(at.x / p);
    }
  }
  Vector AtZPole()
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)Math.cos(theta);
    double sx = (double)Math.sin(theta);
    //
    return new Vector(p * sy * cx, p * sy * sx, p * cy, p);
  }
  Vector AtYPole()
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)Math.cos(theta);
    double sx = (double)Math.sin(theta);
    //
    return new Vector(p * sx * cy, p * sy, p * cx * cy, p);
  }
  Vector AtXPole()
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)Math.cos(theta);
    double sx = (double)Math.sin(theta);
    //
    return new Vector(p * cy, p * sx * sy, p * cx * cy, p);
  }
  Normal ToYPole()
  {

    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)Math.cos(theta);
    double sx = (double)Math.sin(theta);
    Normal ret = new Normal(sx * cy, sy, cx * cy);
    ret.Normalize();
    return ret;
  }
  Normal ToXPole()
  {

    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)Math.cos(theta);
    double sx = (double)Math.sin(theta);
    Normal ret = new Normal(cy, sx * sy, cx * cy);
    ret.Normalize();
    return ret;
  }
  Normal ToZPole()
  {
    double sp = Math.sin(phi);
    //
    Normal n = new Normal(sp * Math.cos(theta), sp * Math.sin(theta), Math.cos(phi));
    n.Normalize();
    return n;
  }
  Normal SphericalDirectionZPole(double sintheta, double costheta, double phi)
  {
    Normal ret = new Normal(sintheta * Math.cos(phi), 
      sintheta * Math.sin(phi), 
      costheta);
    ret.Normalize();
    return ret;
  }
  Normal SphericalDirectionYPole(double sintheta, double costheta, double phi)
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)costheta;
    double sx = (double)sintheta;
    Normal ret = new Normal(sx * cy, sy, cx * cy);
    ret.Normalize();
    return ret;
  }
  Normal SphericalDirectionXPole(double sintheta, double costheta, double phi)
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)costheta;
    double sx = (double)sintheta;
    Normal ret = new Normal(cy, sx * sy, cx * cy);
    ret.Normalize();
    return ret;
  }
  Vector SphericalDirectionZPole(double sintheta, double costheta, double phi, Vector x, Vector y, Vector z)
  {
    Vector xs=x.Scaled((double)(sintheta * Math.cos(phi)));
    Vector ys=y.Scaled((double)(sintheta * Math.sin(phi)));
    Vector zs=z.Scaled(costheta);
    return xs.Add(xs, xs.Add(ys, zs)).AsVector();
  }
  Vector SphericalDirectionYPole(double sintheta, double costheta, double phi, Vector x, Vector y, Vector z)
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)costheta;
    double sx = (double)sintheta;
    //return x * sx * cy + y * sy + z * cx * cy;
    Vector xs=x.Scaled((double)(sx* cy));
    Vector ys=y.Scaled((double)(sy));
    Vector zs=z.Scaled(cx*cy);
    return xs.Add(xs, xs.Add(ys, zs)).AsVector();
  }
  Vector SphericalDirectionXPole(double sintheta, double costheta, double phi, Vector x, Vector y, Vector z)
  {
    double sy = (double)Math.sin(phi);
    double cy = (double)Math.cos(phi);
    double cx = (double)costheta;
    double sx = (double)sintheta;
    //return x * cy + y * sx * sy + z * cx * cy;
    Vector xs=x.Scaled((double)(cy));
    Vector ys=y.Scaled((double)(sx*sy));
    Vector zs=z.Scaled(cx*cy);
    return xs.Add(xs, xs.Add(ys, zs)).AsVector();
  }
  double SphericalThetaZPole(Point v)
  {
    return (double)(Math.acos(Clamp(v.z, -1.0f, 1.0f)));
  }
  double SphericalThetaXPole(Point v)
  {
    return (double)(Math.acos(Clamp(v.x, -1.0f, 1.0f)));
  }
  double SphericalThetaYPole(Point v)
  {
    return (double)(Math.acos(Clamp(v.y, -1.0f, 1.0f)));
  }
  double SphericalPhiZPole(Point v)
  {
    double p = (double)Math.atan2(v.y, v.x);
    return (p < 0.0f) ? p + PI2 : p;
  }
  double SphericalPhiXPole(Point v)
  {
    double p = (double)Math.atan2(v.y, v.z);
    return (p < 0.0f) ? p + PI2 : p;
  }
  double SphericalPhiYPole(Point v)
  {
    double p = (double)Math.atan2(v.z, v.x);
    return (p < 0.0f) ? p + PI2 : p;
  }
}
