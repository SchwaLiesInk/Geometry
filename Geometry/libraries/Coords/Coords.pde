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

