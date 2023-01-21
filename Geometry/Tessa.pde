class Vertex {
  SphericalCoord sv;
  Vector v;
  Vector p;
  Normal n;
  double h;
  int index;
  int adjacent;
  Vertex adj[]=new Vertex[6];
  Tessa parent=null;
  Vertex(Vector vec, Normal norm, Vector colour, double ht) {
    v=vec;
    n=norm;
    p=colour;
    h=ht;
    sv=new SphericalCoord(v, AxisY);
  } 
  Vector MidPoint(Vertex v2) {
    return new Vector((v.x+v2.v.x)*0.5, (v.y+v2.v.y)*0.5, (v.z+v2.v.z)*0.5);
  }
  boolean Equal(Vertex test) {
    return Math.abs(test.v.x-v.x)<E && Math.abs(test.v.y-v.y)<E && Math.abs(test.v.z-v.z)<E;
  }
}
final byte CENTRAL=0;
final byte NORTH=1;
final byte EAST=2;
final byte WEST=3;
abstract interface MakeTessa{
abstract Tessa CreateTessa(Vertex v1, Vertex v2, Vertex v3);
}
class Tessa implements  MakeTessa{

  Tessa parent;
  Vertex vert[]=new Vertex[3];
  Cartesian center;
  Tessa tess[]=null;
  Tessa(Vertex v1, Vertex v2, Vertex v3) {
    vert[0]=v1;
    vert[1]=v2;
    vert[2]=v3;
    center=v1.v.Add(v1.v, v2.v);
    center.Scale(0.5);
    center=center.Add(center, v3.v);
    center.Scale(0.5);
    tess=new Tessa[4];
  }
  Tessa CreateTessa(Vertex v1, Vertex v2, Vertex v3){
    
    return new Tessa( v1, v2, v3);
  }
  
  void GetLevel(ArrayList<Tessa> list, int level) {
    level--;
    if (level<=0 || tess[0].tess[0]==null) {
      list.add(tess[0]);
      list.add(tess[1]);
      list.add(tess[2]);
      list.add(tess[3]);
    } else {
      if (tess[0].tess[0]!=null) {
        tess[0].GetLevel(list, level);
        tess[1].GetLevel(list, level);
        tess[2].GetLevel(list, level);
        tess[3].GetLevel(list, level);
      }
    }
  }
  void Guide(Tessa t) {
    if (tess[0]!=null) {

      tess[1].vert[0].parent=tess[1];
      tess[1].vert[1].parent=tess[1];
      tess[1].vert[2].parent=tess[1];
      //
      tess[2].vert[0].parent=tess[2];
      tess[2].vert[1].parent=tess[2];
      tess[2].vert[2].parent=tess[2];
      //
      tess[3].vert[0].parent=tess[3];
      tess[3].vert[1].parent=tess[3];
      tess[3].vert[2].parent=tess[3];
      //
      tess[0].vert[0].parent=tess[0];
      tess[0].vert[1].parent=tess[0];
      tess[0].vert[2].parent=tess[0];
      //
      tess[0].parent=t;
      tess[1].parent=t;
      tess[2].parent=t;
      tess[3].parent=t;
      tess[3].Guide(tess[3]);
      tess[2].Guide(tess[2]);
      tess[1].Guide(tess[1]);
      tess[0].Guide(tess[0]);
    }
  }
  void Transfer(int level) {
    level--;
    if (level<=0) {
      if (tess[0]!=null) {
        Move(tess[0]);
        Move(tess[1]);
        Move(tess[2]);
        Move(tess[3]);
      }
    } else {
      if (tess[0]!=null) {
        tess[0].Transfer(level);
        tess[1].Transfer(level);
        tess[2].Transfer(level);
        tess[3].Transfer(level);
      }
    }
  }
  void Move(Tessa t) {
    if (t.parent!=null) {
      for (int i=0; i<4; i++) {
        t.parent.tess[i].vert[0].h+=t.vert[0].h;
        t.parent.tess[i].vert[0].h*=0.5;
        t.parent.tess[i].vert[1].h+=t.vert[1].h;
        t.parent.tess[i].vert[1].h*=0.5;
        t.parent.tess[i].vert[2].h+=t.vert[2].h;
        t.parent.tess[i].vert[2].h*=0.5;
      }
    }
  }
  void Flatten(Tessa t) {
    if (t.parent!=null) {
      for (int j=0; j<4; j++) {
        int i=arcRandomInt(4);
        if ((t.parent.tess[i].vert[0].h>0 &&t.parent.tess[i].vert[0].h>0)) {
          t.parent.tess[i].vert[0].h+=t.vert[0].h/t.parent.tess[i].vert[0].h;
          t.parent.tess[i].vert[0].h*=0.5;
        } else if ((t.parent.tess[i].vert[0].h<0 &&t.parent.tess[i].vert[0].h<0)) {
          t.parent.tess[i].vert[0].h-=t.vert[0].h/t.parent.tess[i].vert[0].h;
          t.parent.tess[i].vert[0].h*=0.5;
        }
        if ((t.parent.tess[i].vert[1].h>0 &&t.parent.tess[i].vert[1].h>0)) {
          t.parent.tess[i].vert[1].h+=t.vert[1].h/t.parent.tess[i].vert[1].h;
          t.parent.tess[i].vert[1].h*=0.5;
        } else if ((t.parent.tess[i].vert[1].h<0 &&t.parent.tess[i].vert[1].h<0)) {
          t.parent.tess[i].vert[1].h-=t.vert[1].h/t.parent.tess[i].vert[1].h;
          t.parent.tess[i].vert[1].h*=0.5;
        }
        if ((t.parent.tess[i].vert[2].h>0 &&t.parent.tess[i].vert[2].h>0)) {
          t.parent.tess[i].vert[2].h+=t.vert[2].h/t.parent.tess[i].vert[2].h;
          t.parent.tess[i].vert[2].h*=0.5;
        } else if ((t.parent.tess[i].vert[2].h<0 &&t.parent.tess[i].vert[2].h<0)) {
          t.parent.tess[i].vert[2].h-=t.vert[2].h/t.parent.tess[i].vert[2].h;
          t.parent.tess[i].vert[2].h*=0.5;
        }
      }
    }
  }
  void Flattener(int level) {
    level--;
    if (level>0) {
      if (tess[0]!=null) {
        Flatten(tess[0]);
        Flatten(tess[1]);
        Flatten(tess[2]);
        Flatten(tess[3]);
        if (tess[0]!=null) {
          if (arcRandomInt(3)==0)tess[0].Flattener(level);
          if (arcRandomInt(3)==0)tess[1].Flattener(level);
          if (arcRandomInt(3)==0)tess[2].Flattener(level);
          if (arcRandomInt(3)==0)tess[3].Flattener(level);
        }
      }
    }
  }
  int Equity(Tessa test) {

    //double scale=mass*temp;
    int n=0;
    if (test!=null)
      for (int j=0; j<3; j++) {
        if (test.vert[0].Equal(vert[j])) {
          vert[j]=test.vert[0];
          n++;
        }
        if (test.vert[1].Equal(vert[j])) {
          vert[j]=test.vert[1];
          n++;
        }
        if (test.vert[2].Equal(vert[j])) {
          vert[j]=test.vert[2];
          n++;
        }
      }
    return n;
  }
  void Equality(Tessa test) {
    if (test!=null) {
      if (test.tess[0]!=null && tess[0]!=null) { 

        for (int i=0; i<4; i++) {
          tess[i].Equity(test);
          tess[i].Equality(test);
          for (int j=i+1; j<4; j++) {
            tess[i].Equity(test.tess[j]);
            tess[i].Equality(test.tess[j]);
          }
        }
      }
    }
  }
  void Add(Tessa test) {

    if (tess[0].tess[0]==null) {
      for (int i=0; i<4; i++) {
        for (int j=0; j<3; j++) {
          tess[i].vert[j].h+=test.vert[0].h*Noise(1+(float)test.vert[0].h);
          tess[i].vert[j].p=tess[i].vert[0].p.Add(tess[i].vert[j].p, test.vert[j].p).AsVector();
          tess[i].vert[j].p.Scale(0.5);
        }
      }
    } else {
      tess[0].Add(test.tess[0]);
      tess[1].Add(test.tess[1]);
      tess[2].Add(test.tess[2]);
      tess[3].Add(test.tess[3]);
    }
  }

  void Scale(Tessa test, double scale) {

    if (tess[0].tess[0]==null) {
      for (int i=0; i<4; i++) {
        for (int j=0; j<3; j++) {
          tess[i].vert[j].h*=scale;
        }
      }
    } else {
      tess[0].Scale(test.tess[0], scale);
      tess[1].Scale(test.tess[1], scale);
      tess[2].Scale(test.tess[2], scale);
      tess[3].Scale(test.tess[3], scale);
    }
  }
  void Merge(Tessa t) {

    if (t!=null && tess[0]!=null&& tess[0].tess[0]==null) {
      for (int n=0; n<4; n++) {
        int i=arcRandomInt(4);  

        tess[i].vert[0].h=tess[i].vert[0].h*617+t.vert[2].h;
        tess[i].vert[0].h*=0.00161812297734627831715210355987;
        tess[i].vert[1].h=tess[i].vert[1].h*617+t.vert[0].h;
        tess[i].vert[1].h*=0.00161812297734627831715210355987;
        tess[i].vert[2].h=tess[i].vert[2].h*617+t.vert[1].h;
        tess[i].vert[2].h*=0.00161812297734627831715210355987;
      }
    } else if (t!=null){
      //for (int i=0; i<4; i++) {
      tess[arcRandomInt(4)].Merge(t.tess[arcRandomInt(4)]);
      tess[arcRandomInt(4)].Merge(t.tess[arcRandomInt(4)]);
      tess[arcRandomInt(4)].Merge(t.tess[arcRandomInt(4)]);
      tess[arcRandomInt(4)].Merge(t.tess[arcRandomInt(4)]);
      //}
    }
  }
  void Blur(int level, int to, Tessa t) {
    level--;
    if (level<=to) {
      if (tess[0]!=null  ) {
        for (int i=0; i<4; i++) {
          tess[i].vert[0].h+=t.vert[0].h*t.vert[0].h+t.vert[1].h;
          tess[i].vert[1].h+=t.vert[1].h*t.vert[1].h+t.vert[2].h;
          tess[i].vert[2].h+=t.vert[2].h*t.vert[2].h+t.vert[0].h;

          tess[i].vert[0].h*=0.125;
          tess[i].vert[1].h*=0.125;
          tess[i].vert[2].h*=0.125;
        }
      }
    } else {
      Blur(level, to, t);
    }
  } 
  void ContourLevel() {
    //double ht[]=new double[3];

    for (int i=1; i<4; i++) {
      double total=tess[i].vert[0].h+tess[i].vert[1].h+tess[i].vert[2].h;
      for (int j=0; j<3; j++) {

        tess[i].vert[j].h+=total;//+ht[j];//+(arcRandomAverage(land/water)-arcRandomAverage(water/land))*Math.PI*land;
        tess[i].vert[j].h*=0.25;
      }
    }
  }
}
class PlanetTessa extends Tessa implements MakeTessa{

   
  PlanetTessa(Vertex v1, Vertex v2, Vertex v3) {
    super(v1,v2,v3);
  }
  Tessa CreateTessa(Vertex v1, Vertex v2, Vertex v3){
    
    return (Tessa)new PlanetTessa( v1, v2, v3);
  }
  void SubDivide(Vertex v1, Vertex v2, Vertex v3, GeoPlanet planet, int level) {
    Vector p12=new Vector(v1.MidPoint(v2));
    Vector p23=new Vector(v2.MidPoint(v3));
    Vector p31=new Vector(v3.MidPoint(v1));
    double h1=(v1.h+v2.h)*Noise((float)v1.h, (float)v2.h)*planet.weather+Noise((float)v3.h)*planet.land-Noise((float)v2.h)*planet.water;
    double h2=(v2.h+v3.h)*Noise((float)v2.h, (float)v3.h)*planet.weather+Noise((float)v1.h)*planet.land-Noise((float)v3.h)*planet.water;
    double h3=(v1.h+v3.h)*Noise((float)v1.h, (float)v3.h)*planet.weather+Noise((float)v2.h)*planet.land-Noise((float)v1.h)*planet.water;
    Point o=new Point(1);
    Cartesian v12=o.Sub(p12, o);
    Cartesian v23=o.Sub(p23, o);
    Cartesian v31=o.Sub(p31, o);
    //
    v12.Normalize();
    v23.Normalize();
    v31.Normalize();
    //
    p12=v12.Add(v12, o).AsVector();
    p23=v23.Add(v23, o).AsVector();
    p31=v31.Add(v31, o).AsVector();
    Vector c12=v1.p.Add(v1.p, v2.p).Scaled(0.5).AsVector();
    Vector c23=v2.p.Add(v2.p, v3.p).Scaled(0.5).AsVector();
    Vector c31=v3.p.Add(v3.p, v1.p).Scaled(0.5).AsVector();
    //
    /*int fracLevel=level;
     for (int j=0; j<fracLevel; j++) {
     c12=c12.Add(c12, new Vector(arcRandomAverage(planet.wind), arcRandomAverage(planet.land), arcRandomAverage(planet.water))).AsVector();
     c23=c23.Add(c23, new Vector(arcRandomAverage(planet.wind), arcRandomAverage(planet.land), arcRandomAverage(planet.water))).AsVector();
     c31=c31.Add(c31, new Vector(arcRandomAverage(planet.wind), arcRandomAverage(planet.land), arcRandomAverage(planet.water))).AsVector();
     //
     c12=c12.Add(c12, NoiseVector(p12, planet.mass)).AsVector();
     c23=c23.Add(c23, NoiseVector(p23, planet.mass)).AsVector();
     c31=c31.Add(c31, NoiseVector(p31, planet.mass)).AsVector();
     int rndLevel=level+arcRandomInt(6);
     
     double t=noise((float)center.x, (float)center.z);
     for (int k=1; k<rndLevel; k++) {
     Vector r1=c12.Add(c23, c31).AsVector();
     Vector r2=c23.Add(c12, c31).AsVector();
     Vector r3=c31.Add(c23, c12).AsVector();
     double t1=noise((float)r1.y, (float)r1.z)*planet.soil;
     double t2=noise((float)r2.y, (float)r2.z)*planet.soil;
     double t3=noise((float)r3.y, (float)r3.z)*planet.soil;
     r1.Scale((AverageNoise(c12.x, c12.y, c12.z, 1)+T+Math.abs(p12.y*t1))/(c12.Length()));
     r2.Scale((AverageNoise(c23.x, c23.y, c23.z, 1)+T+Math.abs(p23.y*t2))/(c23.Length()));
     r3.Scale((AverageNoise(c31.x, c31.y, c31.z, 1)+T+Math.abs(p31.y*t3))/(c31.Length()));
     c12=c12.Add(c12, r1).AsVector();
     c23=c23.Add(c23, r2).AsVector();
     c31=c31.Add(c31, r3).AsVector();
     }
     }
     
     c12.Scale(0.25/fracLevel);
     c23.Scale(0.25/fracLevel);
     c31.Scale(0.25/fracLevel);*/

    Clamp(c12);
    Clamp(c23);
    Clamp(c31);
    //
    Vertex vx12=new Vertex(p12, p12.AsNormal(), c12, h1);
    Vertex vx23=new Vertex(p23, p23.AsNormal(), c23, h2);
    Vertex vx31=new Vertex(p31, p31.AsNormal(), c31, h3);
    tess[0]=CreateTessa(vx12, vx23, vx31);
    tess[1]=CreateTessa(v1, vx12, vx31);
    tess[2]=CreateTessa(v2, vx23, vx12);
    tess[3]=CreateTessa(v3, vx31, vx23);
    level--;
    if (level>0) {
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
      t0.SubDivide(tess[0].vert[0], tess[0].vert[1], tess[0].vert[2], planet, level);
      t1.SubDivide(tess[1].vert[0], tess[1].vert[1], tess[1].vert[2], planet, level);
      t2.SubDivide(tess[2].vert[0], tess[2].vert[1], tess[2].vert[2], planet, level);
      t3.SubDivide(tess[3].vert[0], tess[3].vert[1], tess[3].vert[2], planet, level);
    }
  }
  /*void SubDivide(Vertex v1, Vertex v2, Vertex v3, GeoPlanet planet) {
   Vector p12=new Vector(v1.MidPoint(v2));
   Vector p23=new Vector(v2.MidPoint(v3));
   Vector p31=new Vector(v3.MidPoint(v1));
   double h1=(v1.h+v2.h)*noise((float)v1.h, (float)v2.h)*planet.weather+noise((float)v3.h)*planet.land-noise((float)v2.h)*planet.water;
   double h2=(v2.h+v3.h)*noise((float)v2.h, (float)v3.h)*planet.weather+noise((float)v1.h)*planet.land-noise((float)v3.h)*planet.water;
   double h3=(v1.h+v3.h)*noise((float)v1.h, (float)v3.h)*planet.weather+noise((float)v2.h)*planet.land-noise((float)v1.h)*planet.water;
   Point o=new Point(1);
   Cartesian v12=o.Sub(p12, o);
   Cartesian v23=o.Sub(p23, o);
   Cartesian v31=o.Sub(p31, o);
   //
   v12.Normalize();
   v23.Normalize();
   v31.Normalize();
   //
   p12=v12.Add(v12, o).AsVector();
   p23=v23.Add(v23, o).AsVector();
   p31=v31.Add(v31, o).AsVector();
   Vector c12=v1.p.Add(v1.p, v2.p).Scaled(0.5).AsVector();
   Vector c23=v1.p.Add(v2.p, v3.p).Scaled(0.5).AsVector();
   Vector c31=v1.p.Add(v3.p, v1.p).Scaled(0.5).AsVector();
   //
   Vertex vx12=new Vertex(p12, p12.AsNormal(), c12, h1);
   Vertex vx23=new Vertex(p23, p23.AsNormal(), c23, h2);
   Vertex vx31=new Vertex(p31, p31.AsNormal(), c31, h3);
   tess[0]=new Tessa(vx12, vx23, vx31);
   tess[1]=new Tessa(v1, vx12, vx31);
   tess[2]=new Tessa(v2, vx23, vx12);
   tess[3]=new Tessa(v3, vx31, vx23);
   }*/
  void ColorizeLevel(GeoPlanet planet, Vector stain) {
    Vector c=new Vector(AverageNoise(center.x* planet.atmos, center.y* planet.atmos, center.z* planet.atmos, planet.clouds), 
      AverageNoise(center.x*planet.water, center.y*planet.water, center.z*planet.water, planet.rivers), 
      AverageNoise(planet.land*center.z, planet.land*center.y, planet.land*center.x, planet.rain));

    for (int i=0; i<4; i++) {
      for (int j=0; j<3; j++) {

        double ht=tess[i].vert[j].h/planet.mass;
        //Clamp(ht,-water,rain);
        double lum=(1+Math.abs(tess[i].vert[j].v.y))*planet.temp;
        if (tess[i].vert[j].p.y>tess[i].vert[j].p.z) {
          tess[i].vert[j].p.z*=lum*(planet.rain/planet.temp);
        }

        tess[i].vert[j].p.z+=Math.abs(ht*planet.water);
        if (tess[i].vert[j].p.z>tess[i].vert[j].p.y) {
          tess[i].vert[j].p.y*=lum*(planet.rivers/planet.temp);
        }
        if (tess[i].vert[j].h>0) {
          tess[i].vert[j].p.y+=ht*planet.water;
        }
        if (Math.abs(tess[i].vert[j].v.y)>tess[i].vert[j].p.x) {
          tess[i].vert[j].p.x*=lum*(planet.temp/planet.clouds);
        }
        tess[i].vert[j].p.x+=ht/planet.atmos;

        if (tess[i].vert[j].p.x<0) {
          tess[i].vert[j].p.x=0;
        }
        c=c.Add(c, new Vector(AverageNoise(1+center.z, 1+center.x, 1+center.y, planet.wind), 
          AverageNoise(1+center.x, 1+center.y, 1+center.z, planet.land), 
          AverageNoise(1+center.z, 1+center.y, 1+center.x, planet.water))).AsVector();

        //c=c.Add(c, new Vector(arcRandomAverage(wind), arcRandomAverage(land), arcRandomAverage(rain))).AsVector();
        c.Scale(0.5);
        tess[i].vert[j].p=c.Add(c, tess[i].vert[j].p).AsVector();
        tess[i].vert[j].p.Scale(0.5);
        //if (Math.abs(tess[i].vert[j].v.y)>tess[i].vert[j].p.z) {
        lum=(2.73-planet.temp)*Math.abs(tess[i].vert[j].v.y)*Math.abs(tess[i].vert[j].v.y);
        tess[i].vert[j].p.x+=lum*(planet.clouds/planet.temp)/planet.size;//+ht;
        tess[i].vert[j].p.y+=lum*(planet.rivers/planet.temp)*planet.size;//+ht;
        tess[i].vert[j].p.z+=lum*(planet.temp/planet.rain)*planet.size;//+Math.abs(ht);
        //}
      }
      c=c.Add(c, tess[i].vert[0].p).AsVector();
      c=c.Add(c, tess[i].vert[1].p).AsVector();
      c=c.Add(c, tess[i].vert[2].p).AsVector();
      c.Scale(0.25*planet.temp);
      Vector c1=tess[i].vert[0].p=c.Add(c, tess[i].vert[0].p).AsVector();
      tess[i].vert[0].p.Scale(Math.abs(tess[i].vert[0].h));
      Vector c2=tess[i].vert[1].p=c.Add(c, tess[i].vert[1].p).AsVector();
      tess[i].vert[1].p.Scale(Math.abs(tess[i].vert[1].h));
      Vector c3=tess[i].vert[2].p=c.Add(c, tess[i].vert[2].p).AsVector();
      tess[i].vert[2].p.Scale(Math.abs(tess[i].vert[2].h));
      //
      Vector r1=c1.Add(c2, c3).AsVector();
      Vector r2=c2.Add(c1, c3).AsVector();
      Vector r3=c3.Add(c2, c1).AsVector();
      double t1=Noise((float)r1.y, (float)r1.z)*(Math.PI-planet.temp+Noise((float)r1.x, (float)r1.z)*0.1);
      double t2=Noise((float)r2.y, (float)r2.z)*(Math.PI-planet.temp+Noise((float)r2.x, (float)r2.z)*0.1);
      double t3=Noise((float)r3.y, (float)r3.z)*(Math.PI-planet.temp+Noise((float)r3.x, (float)r3.z)*0.1);
      r1.Scale((AverageNoise(c1.x, c1.y, c1.z, 1)+T+Math.abs(tess[i].vert[0].v.y*t1))/(1+c1.Length()));
      r2.Scale((AverageNoise(c2.x, c2.y, c2.z, 1)+T+Math.abs(tess[i].vert[1].v.y*t2))/(1+c2.Length()));
      r3.Scale((AverageNoise(c3.x, c3.y, c3.z, 1)+T+Math.abs(tess[i].vert[2].v.y*t3))/(1+c3.Length()));
      tess[i].vert[0].p=c1.Add(c1, r1).AsVector();
      tess[i].vert[1].p=c2.Add(c2, r2).AsVector();
      tess[i].vert[2].p=c3.Add(c3, r3).AsVector();
      tess[i].vert[0].p.Scale(T+0.25*planet.temp);
      tess[i].vert[1].p.Scale(T+0.25*planet.temp);
      tess[i].vert[2].p.Scale(T+0.25*planet.temp);
      //
      double cy=(T+Math.abs(0.5-Math.abs(tess[i].center.y)))*(planet.temp);
      double cx=Math.abs(T+1+Math.abs(tess[i].center.y));
      tess[i].vert[0].p=tess[i].vert[0].p.Add(tess[i].vert[0].p, new Vector(tess[i].vert[0].h, tess[i].vert[0].h, tess[i].vert[0].h)).AsVector();
      tess[i].vert[1].p=tess[i].vert[1].p.Add(tess[i].vert[1].p, new Vector(tess[i].vert[1].h, tess[i].vert[1].h, tess[i].vert[1].h)).AsVector();
      tess[i].vert[2].p=tess[i].vert[2].p.Add(tess[i].vert[2].p, new Vector(tess[i].vert[2].h, tess[i].vert[2].h, tess[i].vert[2].h)).AsVector();
      tess[i].vert[0].p.Scale((cx+cy)/planet.mass);
      tess[i].vert[1].p.Scale((cx+cy)/planet.mass);
      tess[i].vert[2].p.Scale((cx+cy)/planet.mass);
      if (tess[i].vert[0].p.z<Epsilon) {
        tess[i].vert[0].p.z+=Math.abs(tess[i].vert[0].h);
        tess[i].vert[0].p.z=1.0-Noise(1+(float)tess[i].vert[0].p.z);
      }
      if (tess[i].vert[1].p.z<Epsilon) {
        tess[i].vert[1].p.z+=Math.abs(tess[i].vert[1].h);
        tess[i].vert[1].p.z=1.0-Noise(1+(float)tess[i].vert[1].p.z);
      }
      if (tess[i].vert[2].p.z<Epsilon) {
        tess[i].vert[2].p.z+=Math.abs(tess[i].vert[2].h);
        tess[i].vert[2].p.z=1.0-Noise(1+(float)tess[i].vert[2].p.z);
      }
      if (tess[i].vert[0].h>planet.atmos*2) {
        tess[i].vert[0].p.x/=tess[i].vert[0].h;
        tess[i].vert[0].p.y/=tess[i].vert[0].h;
        tess[i].vert[0].p.z/=tess[i].vert[0].h;
      }
      if (tess[i].vert[1].h>planet.atmos*2) {
        tess[i].vert[1].p.x/=tess[i].vert[1].h;
        tess[i].vert[1].p.y/=tess[i].vert[1].h;
        tess[i].vert[1].p.z/=tess[i].vert[1].h;
      }
      if (tess[i].vert[2].h>planet.atmos*2) {
        tess[i].vert[2].p.x/=tess[i].vert[2].h;
        tess[i].vert[2].p.y/=tess[i].vert[2].h;
        tess[i].vert[2].p.z/=tess[i].vert[2].h;
      }
      //
      if (tess[i].vert[0].h<planet.atmos && tess[i].vert[0].h>0) {
        tess[i].vert[0].p.y+=planet.rain/planet.mass;
      }
      if (tess[i].vert[1].h<planet.atmos && tess[i].vert[1].h>0) {
        tess[i].vert[1].p.y+=planet.rain/planet.mass;
      }
      if (tess[i].vert[2].h<planet.atmos && tess[i].vert[2].h>0) {
        tess[i].vert[2].p.y+=planet.rain/planet.mass;
      }

      if (tess[i].vert[0].h<planet.atmos && tess[i].vert[0].h>-planet.atmos) {
        tess[i].vert[0].p.y+=planet.rain/planet.mass;
      }
      if (tess[i].vert[1].h<planet.atmos && tess[i].vert[1].h>-planet.atmos) {
        tess[i].vert[1].p.y+=planet.rain/planet.mass;
      }
      if (tess[i].vert[2].h<planet.atmos && tess[i].vert[2].h>-planet.atmos) {
        tess[i].vert[2].p.y+=planet.rain/planet.mass;
      }
      //
      if (tess[i].vert[0].h<-planet.water) {
        tess[i].vert[0].p.y+=planet.water/Math.PI;
      }
      if (tess[i].vert[1].h<-planet.water) {
        tess[i].vert[1].p.y+=planet.water/Math.PI;
      }
      if (tess[i].vert[2].h<-planet.water) {
        tess[i].vert[2].p.y+=planet.water/Math.PI;
      }

      if (tess[i].vert[0].p.x<planet.temp) {
        tess[i].vert[0].p.x+=planet.wind/planet.mass;
      }
      if (tess[i].vert[1].p.x<planet.temp) {
        tess[i].vert[1].p.x+=planet.wind/planet.mass;
      }
      if (tess[i].vert[2].p.x<planet.temp) {
        tess[i].vert[2].p.x+=planet.wind/planet.mass;
      }
      c1=NoiseVector(tess[i].vert[0].p, 2);
      c2=NoiseVector(tess[i].vert[1].p, 2);
      c3=NoiseVector(tess[i].vert[2].p, 2);
      c1=c1.Add(c1, stain).AsVector();
      c1.Scale(0.25+(Math.abs(Math.asin(tess[i].vert[0].v.y)))/(planet.temp*Math.PI));
      //if(tess[i].vert[0].h>0)c1.Scale(Math.abs(tess[i].vert[0].h*planet.clouds)/(planet.temp*Math.PI*4));
      tess[i].vert[0].p=c1.Add(c1, tess[i].vert[0].p).AsVector();
      tess[i].vert[0].p.Scale(0.5);
      c2=c2.Add(c2, stain).AsVector();
      c2.Scale(0.25+(Math.abs(Math.asin(tess[i].vert[1].v.y)))/(planet.temp*Math.PI));
      //if(tess[i].vert[1].h>0)c2.Scale(Math.abs(tess[i].vert[1].h*planet.clouds)/(planet.temp*Math.PI*4));
      tess[i].vert[1].p=c2.Add(c2, tess[i].vert[1].p).AsVector();
      tess[i].vert[1].p.Scale(0.5);
      c3=c3.Add(c3, stain).AsVector();
      c3.Scale(0.25+(Math.abs(Math.asin(tess[i].vert[2].v.y)))/(planet.temp*Math.PI));
      //if(tess[i].vert[2].h>0)c3.Scale((tess[i].vert[2].h*planet.clouds)/(planet.temp*Math.PI*4));
      tess[i].vert[2].p=c3.Add(c3, tess[i].vert[2].p).AsVector();
      tess[i].vert[2].p.Scale(0.5);
      if (tess[i].vert[0].h<0) {
        tess[i].vert[0].p.z+=planet.water/(planet.mass/planet.temp);
      } else {
        tess[i].vert[0].p.y+=planet.rain/(planet.mass+tess[i].vert[0].h);
      }
      if (tess[i].vert[1].h<0) {
        tess[i].vert[1].p.z+=planet.water/(planet.mass/planet.temp);
      } else {
        tess[i].vert[1].p.y+=planet.rain/(planet.mass+tess[i].vert[1].h);
      }
      if (tess[i].vert[2].h<0) {
        tess[i].vert[2].p.z+=planet.water/(planet.mass/planet.temp);
      } else {
        tess[i].vert[2].p.y+=planet.rain/(planet.mass+tess[i].vert[2].h);
      }
      Clamp(tess[i].vert[0].p);
      Clamp(tess[i].vert[1].p);
      Clamp(tess[i].vert[2].p);
    }
  }
  void ColorizeLevel(GeoPlanet planet, Vector stain, int level) {

    level--;
    if (level<=0 || tess[0].tess==null|| tess[0].tess[0]==null) {
      ColorizeLevel(planet, stain);
    } else {
      //ColorizeLevel(planet);
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
      t0.ColorizeLevel(planet, stain, level);
      t1.ColorizeLevel(planet, stain, level);
      t2.ColorizeLevel(planet, stain, level);
      t3.ColorizeLevel(planet, stain, level);
    }
  }
 
  void ContourLevel(GeoPlanet planet, int level) {
    level--;
    if (level<=0) {
      if (tess[0]!=null  ) {

        for (int i=0; i<4; i++) {

          double total=tess[i].vert[0].h+tess[i].vert[1].h+tess[i].vert[2].h;
          for (int j=0; j<3; j++) {

            tess[i].vert[j].h+=total;//+ht[j];//+(arcRandomAverage(land/water)-arcRandomAverage(water/land))*Math.PI*land;

            tess[i].vert[j].h*=0.25;
            if (tess[i].vert[j].h>planet.atmos*Math.PI*2) {
              tess[i].vert[j].h/=Math.PI;
            }
            if (tess[i].vert[j].h>planet.atmos*Math.PI) {
              tess[i].vert[j].h-=Noise((float)(tess[i].vert[j].h*planet.wind), (float)(tess[i].vert[j].h*planet.rain))*planet.clouds*planet.size;
            }
          }
        }
      }
    } else {
      ContourLevel(planet, level);
    }
  }

  void RiverLevel(int level, int to) {
    level--;
    if (level<=to) {
      if (tess[0]!=null) {
        int highi=0;
        int highj=0;
        double highest=-MAX_MAP;
        int lowi=0;
        int lowj=0;
        double lowest=MAX_MAP;
        for (int i=0; i<4; i++) {

          for (int j=0; j<3; j++) {

            if ( tess[i].vert[j].h<lowest) {
              lowest=tess[i].vert[j].h;
              lowi=i;
              lowj=j;
            }
            if ( tess[i].vert[j].h>highest) {
              highest=tess[i].vert[j].h;
              highi=i;
              highj=j;
            }
          }
        }
        if ((highest>0 && lowest>0) || (highest<0 && lowest<0)) {
          tess[highi].vert[highj].h+=tess[lowi].vert[lowj].h/tess[highi].vert[highj].h;
          tess[highi].vert[highj].h*=0.5;
          tess[lowi].vert[lowj].h=tess[highi].vert[highj].h;
          if ((highest>0 && lowest>0)) {
            double w1=Math.abs(tess[lowi].vert[1].h+tess[lowi].vert[2].h)*0.5;
            double w2=Math.abs(tess[lowi].vert[0].h+tess[lowi].vert[2].h)*0.5;
            double w3=Math.abs(tess[lowi].vert[1].h+tess[lowi].vert[0].h)*0.5;
            for (int i=0; i<4; i++) {
              if (tess[i].vert[0].h>w1) {
                tess[i].vert[0].h+=w1;
                tess[i].vert[0].h*=0.5;
              }
              if (tess[i].vert[1].h>w2) {
                tess[i].vert[1].h+=w2;
                tess[i].vert[1].h*=0.5;
              }
              if (tess[i].vert[2].h>w3) {
                tess[i].vert[2].h+=w3;
                tess[i].vert[2].h*=0.5;
              }
            }
            if (this.parent!=null) {
              if (this.parent.vert[0].h>w1) {
                this.parent.vert[0].h+=w1;
                this.parent.vert[0].h*=0.5;
              }
              if (this.parent.vert[1].h>w2) {
                this.parent.vert[1].h+=w2;
                this.parent.vert[1].h*=0.5;
              }
              if (this.parent.vert[2].h>w3) {
                this.parent.vert[2].h+=w3;
                this.parent.vert[2].h*=0.5;
              }
            }
          }
        }
      } else

      {
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
        t0.RiverLevel(level, to);
        t1.RiverLevel(level, to);
        t2.RiverLevel(level, to);
        t3.RiverLevel(level, to);
      }
    }
  }

  void Blend(GeoPlanet planet, Tessa t) {
    if (t.parent!=null) {
      double me=1.0/planet.mass;
      for (int i=0; i<4; i++) {
        t.parent.tess[i].vert[0].h+=-t.vert[0].h*me;
        t.parent.tess[i].vert[1].h+=-t.vert[1].h*me;
        t.parent.tess[i].vert[2].h+=-t.vert[2].h*me;
      }
    }
  }
  void Blender(GeoPlanet planet, int level) {
    level--;
    if (level<=0) {
      if (tess[0]!=null) {
        Blend(planet, tess[0]);
        Blend(planet, tess[1]);
        Blend(planet, tess[2]);
        Blend(planet, tess[3]);
      }
    } else {
      if (tess[0]!=null) {
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
        t0.Blender(planet, level);
        t1.Blender(planet, level);
        t2.Blender(planet, level);
        t3.Blender(planet, level);
      }
    }
  }

  void LandFall(Tessa t, double imass) {
    if (parent!=null) {
      double me= imass;
      if ((parent.vert[0].h>0 &&t.vert[0].h>0)) {
        double d=(t.vert[0].h-parent.vert[0].h);
        d*=me;

        t.vert[0].h+=d;
      }
      if ((parent.vert[1].h>0 &&t.vert[1].h>0)) {
        double d=(t.vert[1].h-parent.vert[1].h);
        d*=me;
        t.vert[1].h+=d;
      }
      if ((parent.vert[2].h>0 &&t.vert[2].h>0)) {
        double d=(t.vert[2].h-parent.vert[2].h);
        d*=me;

        t.vert[2].h+=d;
      }
      //}
    }
  }
  void LandFalls(int level, double imass, int depth) {
    level--;
    if (level>=0) {
      if (tess[0]!=null) {
        LandFall(tess[0], imass);
        LandFall(tess[1], imass);
        LandFall(tess[2], imass);
        LandFall(tess[3], imass);
        if (tess[0]!=null) {
          for (int i=0; i<depth-level; i++) {  
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
            if (arcRandomInt(2)==0)t0.LandFalls(level, imass, depth);
            if (arcRandomInt(2)==0)t1.LandFalls(level, imass, depth);
            if (arcRandomInt(2)==0)t2.LandFalls(level, imass, depth);
            if (arcRandomInt(2)==0)t3.LandFalls(level, imass, depth);
          }
        }
      }
    }
  }
  void River(GeoPlanet planet, Tessa t) {
    if (parent!=null) {

      if ((parent.vert[0].h>0 &&t.vert[0].h>0)) {
        double d=Math.abs(t.vert[0].h-parent.vert[0].h);
        d*=planet.mass;

        t.vert[0].h-=d*(planet.rain+planet.rivers);
        t.vert[0].h*=0.5;
      }
      if ((parent.vert[1].h>0 &&t.vert[1].h>0)) {
        double d=Math.abs(t.vert[1].h-parent.vert[1].h);
        d*=planet.mass;

        t.vert[1].h-=d*(planet.rain+planet.rivers);
        t.vert[1].h*=0.5;
      }
      if ((parent.vert[2].h>0 &&t.vert[2].h>0)) {
        double d=Math.abs(t.vert[2].h-parent.vert[2].h);
        d*=planet.mass;


        t.vert[2].h-=d*(planet.rain+planet.rivers);
        t.vert[2].h*=0.5;
      }
      //}
    }
  }
  void Riverer(GeoPlanet planet, int level) {
    level--;
    if (tess[0]!=null) {
      double av0=(vert[0].h+vert[1].h+vert[2].h);
      if (av0>Epsilon) {
        double av1=(tess[0].vert[0].h+tess[0].vert[1].h+tess[0].vert[2].h);
        double av2=(tess[1].vert[0].h+tess[1].vert[1].h+tess[1].vert[2].h);
        double av3=(tess[2].vert[0].h+tess[2].vert[1].h+tess[2].vert[2].h);
        double av4=(tess[3].vert[0].h+tess[3].vert[1].h+tess[3].vert[2].h);
        if (av1<av0)River(planet, tess[0]);
        if (av2<av0)River(planet, tess[1]);
        if (av3<av0)River(planet, tess[2]);
        if (av4<av0)River(planet, tess[3]);
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
        if (av1<av0)t0.Riverer(planet, level);
        if (av2<av0)t1.Riverer(planet, level);
        if (av3<av0)t2.Riverer(planet, level);
        if (av4<av0)t3.Riverer(planet, level);
      }
    }
  }
  void Seed(int level, int from, int to, double seed) {
    level--;
    if (level<from && level>to) {
      if (tess[0]!=null  ) {
        double av=vert[0].h+vert[1].h+vert[2].h
          +tess[0].vert[0].h+tess[0].vert[1].h+tess[0].vert[2].h
          +tess[1].vert[0].h+tess[1].vert[1].h+tess[1].vert[2].h
          +tess[2].vert[0].h+tess[2].vert[1].h+tess[2].vert[2].h
          +tess[3].vert[0].h+tess[3].vert[1].h+tess[3].vert[2].h;
        if ((av>0 && seed>0) || (av<0 && seed<0)) {
          seed=-seed*Math.abs(av/15);
        }
        for (int i=0; i<4; i++) {

          double y=(1+Math.abs(tess[i].center.y));
          tess[i].vert[0].h+=seed*y*Noise((float)vert[1].h, (float)vert[2].h);
          tess[i].vert[1].h+=seed*y*Noise((float)vert[0].h, (float)vert[2].h);
          tess[i].vert[2].h+=seed*y*Noise((float)vert[1].h, (float)vert[0].h);

          tess[i].vert[0].h*=0.5;
          tess[i].vert[1].h*=0.5;
          tess[i].vert[2].h*=0.5;
        }
      }
    } else {
      if (level>=0 && tess[0]!=null  ) {

        double av=vert[0].h+vert[1].h+vert[2].h
          +tess[0].vert[0].h+tess[0].vert[1].h+tess[0].vert[2].h
          +tess[1].vert[0].h+tess[1].vert[1].h+tess[1].vert[2].h
          +tess[2].vert[0].h+tess[2].vert[1].h+tess[2].vert[2].h
          +tess[3].vert[0].h+tess[3].vert[1].h+tess[3].vert[2].h;
        if ((av>0 && seed>0) || (av<0 && seed<0)) {
          seed=-seed*Math.abs(av/15);
        }  
      PlanetTessa t0=(PlanetTessa)(tess[0]);
      PlanetTessa t1=(PlanetTessa)(tess[1]);
      PlanetTessa t2=(PlanetTessa)(tess[2]);
      PlanetTessa t3=(PlanetTessa)(tess[3]);
        t0.Seed(arcRandomInt(1+level), from, to, seed); 
        t1.Seed(arcRandomInt(1+level), from, to, seed); 
        t2.Seed(arcRandomInt(1+level), from, to, seed); 
        t3.Seed(arcRandomInt(1+level), from, to, seed);
      }
    }
  }
}
