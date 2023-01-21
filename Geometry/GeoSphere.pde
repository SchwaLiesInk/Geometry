
final double POLYX =0.525731112119133606;
final double  POLYZ =0.850650808352039932;
final double vdata[][]= {

  { 
    -POLYX, 0, POLYZ
  }
  , { 
    POLYX, 0, POLYZ
  }
  , { 
    -POLYX, 0, -POLYZ
  }
  , { 
    POLYX, 0, -POLYZ
  }
  , 
  { 
    0, POLYZ, POLYX
  }
  , { 
    0, POLYZ, -POLYX
  }
  , { 
    0, -POLYZ, POLYX
  }
  , { 
    0, -POLYZ, -POLYX
  }
  , 
  { 
    POLYZ, POLYX, 0
  }
  , { 
    -POLYZ, POLYX, 0
  }
  , { 
    POLYZ, -POLYX, 0
  }
  , { 
    -POLYZ, -POLYX, 0
  }
};
final byte indexData[][]= {
  { 
    0, 4, 1
  }
  , { 
    0, 9, 4
  }
  , { 
    9, 5, 4
  }
  , { 
    4, 5, 8
  }
  , { 
    4, 8, 1
  }
  , 
  { 
    8, 10, 1
  }
  , { 
    8, 3, 10
  }
  , { 
    5, 3, 8
  }
  , { 
    5, 2, 3
  }
  , { 
    2, 7, 3
  }
  , 
  { 
    7, 10, 3
  }
  , { 
    7, 6, 10
  }
  , { 
    7, 11, 6
  }
  , { 
    11, 0, 6
  }
  , { 
    0, 1, 6
  }
  , 
  { 
    6, 1, 10
  }
  , { 
    9, 0, 11
  }
  , { 
    9, 11, 2
  }
  , { 
    9, 2, 5
  }
  , { 
    7, 2, 11
  }
};

/*
Vector NoiseVector(Vector v1, double scale) {
 mole2.SingleCalc();
 mole1.SingleCalc();
 S=(float)(arcRandomAverage((mole1.Cos()+mole1.Sin()))*scale);
 S*=S;
 return new Vector(arcRandom()*(((1+v1.x)*S)*(1+mole2.Cos())*S),
 arcRandom()*((1+v1.y)*S* (1+mole2.Sin())*S),
 arcRandom()*((1+v1.z)*S* ((2+mole2.Cos()+mole2.Sin())*S)));
 }*/


void DrawPlanet() {

  int sz= geolist.size();
  beginShape(TRIANGLES);
  for (int i=0; i<sz; i++) {
    Tessa draw=geolist.get(i);
    Vector col1=draw.vert[0].p;
    Vector col2=draw.vert[1].p;
    Vector col3=draw.vert[2].p;
    Normal norm1=draw.vert[0].n; 
    Normal norm2=draw.vert[1].n; 
    Normal norm3=draw.vert[2].n; 
    Vector v1=draw.vert[0].v; 
    Vector v2=draw.vert[1].v; 
    Vector v3=draw.vert[2].v;    
    fill((float)col2.x*255, (float) col2.y*255, (float) col2.z*255);
    normal((float)norm2.x, (float)norm2.y, (float)norm2.z);
    vertex((float)v2.x, (float)v2.y, (float) v2.z);  
    fill((float)col1.x*255, (float) col1.y*255, (float) col1.z*255);
    normal((float)norm1.x, (float)norm1.y, (float)norm1.z);
    vertex((float)v1.x, (float)v1.y, (float)v1.z);
    fill((float)col3.x*255, (float)col3.y*255, (float)col3.z*255);
    normal((float)norm3.x, (float)norm3.y, (float)norm3.z);
    vertex((float)v3.x, (float)v3.y, (float)v3.z);
  }
  endShape();
}
class GeoPlanet {
  StarSystem system;
  StellarBody body;
  double size=1.0;
  double water=0.51;
  double land=1.0-water;
  double temp=0.5;
  double atmos=1.0;
  double mass;
  double rivers;
  double soil;
  double rain;
  double wind;
  double clouds;
  double weather;
  double energy;
  Vector stain;
  GeoPlanet(StarSystem syst) {
    system=syst;
    body=new StellarBody(syst);
    size=(2.0*arcRandom()+(arcRandom()-arcRandom())*arcRandom()*0.5+1)*0.5;
    size=Math.sqrt(size);
    if (size<0.2) {
      size=0.2+Math.abs(size);
    }
    atmos=1.0*arcRandom()+(arcRandom()*size-arcRandom())*arcRandom()*0.5+0.5;
    atmos=Math.sqrt(atmos);
    if (atmos<0.1) {
      atmos=0.1+Math.abs(atmos);
    }
    water=1.0*arcRandom()+((arcRandom()*atmos)-(arcRandom()))*arcRandom()*0.5+0.5;
    water=Math.sqrt(water);
    if (water<0.1) {
      water=0.1+Math.abs(water);
    }
    land=Math.sqrt((size+water)*Math.PI*4)/PI2;
    if (land<0) {
      land=0;
    }
    temp=0.273+3.0*arcRandom()+((arcRandom()*atmos)-(arcRandom()*water))*arcRandom();
    temp=Math.sqrt(temp);
    if (temp<0) {
      temp=0;
    }
    //
    rivers=land*water;
    wind=temp*atmos;
    weather=atmos*size;
    clouds=atmos*water;
    rain=clouds*temp;
    soil=land*rivers;
    mass=(clouds+soil)*Math.PI*4;
    energy=mass*(temp/clouds/size);

    double r=(2.0-size*water)*atmos*temp*arcRandom()+arcRandom()*0.25;
    double g=(land+water)*atmos*arcRandom()+arcRandom()*0.25;
    double b=(water*2)*atmos*arcRandom()+arcRandom()*0.25;
    stain=new Vector(r, g, b);
    println("S"+(int)(size*10)+"W"+(int)(water*10)+"L"+(int)(land*10)+"A"+(int)(atmos*10)+"T"+(int)(temp*10)+"M"+(int)(mass*10));
  
  }
  
  /*GeoPlanet(StellarBody sb) {
    
    size=sb.pr ;//(2.0*arcRandom()+(arcRandom()-arcRandom())*arcRandom()*0.5+1)*0.5;
    //size=Math.sqrt(size);
    //if (size<0.2) {
      //size=0.2+Math.abs(size);
    //}
    double me=sb.m/(Math.PI*2);
    atmos=1.0*arcRandom()+(arcRandom()*size-arcRandom())*arcRandom()*0.5+0.5;
    atmos=Math.sqrt(atmos)*me;
    if (atmos<0.1) {
      atmos=0.1+Math.abs(atmos);
    }
    water=1.0*arcRandom()+((arcRandom()*atmos)-(arcRandom()))*arcRandom()*0.5+0.5;
    water=Math.sqrt(water)*me;
    if (water<0.1) {
      water=0.1+Math.abs(water);
    }
    land=Math.sqrt((size+water)*Math.PI*4)/PI2;
    if (land<0) {
      land=0;
    }
    temp=0.273+3.0*arcRandom()+((arcRandom()*atmos)-(arcRandom()*water))*arcRandom();
    temp=Math.sqrt(temp)*me;
    if (temp<0) {
      temp=0;
    }
    //
    rivers=land*water;
    wind=temp*atmos;
    weather=atmos*size;
    clouds=atmos*water;
    rain=clouds*temp;
    soil=land*rivers;
    mass=sb.m;//(clouds+soil)*Math.PI*4;
    energy=mass*(temp/clouds/size);

    double r=(2.0-size*water)*atmos*temp*arcRandom()+arcRandom()*0.25;
    double g=(land+water)*atmos*arcRandom()+arcRandom()*0.25;
    double b=(water*2)*atmos*arcRandom()+arcRandom()*0.25;
    stain=new Vector(r, g, b);
    println("S"+(int)(size*10)+"W"+(int)(water*10)+"L"+(int)(land*10)+"A"+(int)(atmos*10)+"T"+(int)(temp*10)+"M"+(int)(mass*10));
  }*/
}
class GraphicsOrder {
  ArrayList<Tessa>vertecies=new ArrayList<Tessa>();
}
class GeoSphere {
  GraphicsOrder ge[];
  GeoPlanet planet;
  PlanetTessa triangle[]=new PlanetTessa[20];
  PlanetTessa triangleBuffer[]=new PlanetTessa[20];
  Vector pos=new Vector(0, 0, 0);
  int depth=1;
  GeoSphere(int detailLevel,StarSystem syst) {
    planet=new GeoPlanet(syst);
    if (detailLevel<1) {
      detailLevel=1;
    }
    depth=detailLevel;
    SetUp();
    SetVerts();

    int m=1+(int)(planet.mass);

    FirstSeed();
    SecondSeed(m);
    int w=1+(int)(planet.rain+planet.rivers);
    int age=m*w;
    for (int i=0; i<age; i++) {
      int j=arcRandomInt(12);
      switch(j) {
      case 0:
        { 
          Path(age);
          break;
        }
      case 1:
        { 
          Diffusion(w, 240);
          break;
        }
      case 2:
        { 
          River(); 
          break;
        }
      case 3:
        { 
          Merge(w); 
          break;
        }
      case 4:
        { 
          Rivers(w);
          break;
        }
      case 5:
        { 
          Blur();
          break;
        }
      case 6:
        { 
          Contour();
          break;
        }
      case 7:
        { 
          Transfer();
          break;
        }
      case 8:
        { 
          Flatten();
          break;
        }
      case 9:
        { 
          LandFalls(w);
          break;
        }
      case 10:
        { 
          SecondSeed(w);
          break;
        }
      case 11:
        { 
          SecondSeed(m);
          break;
        }
      }
    }
    Fracture(age>>1);
    /*double r=(2.0-planet.size*planet.water)*planet.atmos*planet.temp;
     double g=(planet.land+planet.water)*planet.atmos;
     double b=(planet.water*2)*planet.atmos;
     Vector c1=new Vector(r,g,b);*/
    Colorize(planet.stain);
  }
  void FirstSeed() {

    for (int j=0; j<120; j++) {
      int i=arcRandomInt(20);
      double y=(1+Math.abs(triangle[i].center.y))*planet.energy*planet.rain*planet.rivers;
      double h1=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*y;
      double h2=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*y;
      double h3=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*y;
      int a1=arcRandomInt(3);
      int a2=arcRandomInt(3);
      int a3=arcRandomInt(3);
      if ((triangle[i].vert[a1].h>0 && h1>0) ||(triangle[i].vert[a1].h<0 && h1<0)) {
        triangle[i].vert[a1].h-=h1;
      } else {
        triangle[i].vert[a1].h+=h1;
      }
      triangle[i].vert[a1].h*=0.5;
      if ((triangle[i].vert[a2].h>0 && h2>0) ||(triangle[i].vert[a2].h<0 && h2<0)) {
        triangle[i].vert[a2].h-=h2;
      } else {
        triangle[i].vert[a2].h+=h2;
      }
      triangle[i].vert[a2].h*=0.5;
      if ((triangle[i].vert[a3].h>0 && h3>0) ||(triangle[i].vert[a3].h<0 && h3<0)) {
        triangle[i].vert[a3].h-=h3;
      } else {
        triangle[i].vert[a3].h+=h3;
      }
      triangle[i].vert[a3].h*=0.5;
    }
  }
  void SecondSeed(int m) {

    for (int n=0; n<m; n++) {
      double seed=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*planet.energy;
      seed+=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*planet.energy;
      seed*=0.5;
      int n8=(n%8);
      //int n0=1+arcRandomInt(n8+1);
      for (int j=0; j<20; j++) {

        double av=triangle[j].vert[0].h+triangle[j].vert[1].h+triangle[j].vert[2].h
          +triangle[j].tess[0].vert[0].h+triangle[j].tess[0].vert[1].h+triangle[j].tess[0].vert[2].h
          +triangle[j].tess[1].vert[0].h+triangle[j].tess[1].vert[1].h+triangle[j].tess[1].vert[2].h
          +triangle[j].tess[2].vert[0].h+triangle[j].tess[2].vert[1].h+triangle[j].tess[2].vert[2].h
          +triangle[j].tess[3].vert[0].h+triangle[j].tess[3].vert[1].h+triangle[j].tess[3].vert[2].h;
        if (av>0) {
          triangle[j].Seed(n8, n8, 0, -Math.abs(seed+av/15)*planet.rain*planet.rivers);
        } else {
          triangle[j].Seed(n8, n8, 0, Math.abs(seed-av/15)*planet.land);
        }
      }
    }
  }
  void SetUp() {

    Vertex vn[][]=new Vertex[20][3];
    Matrix mat0=new Matrix();
    mat0.RotateX(1.0471975511965977461542144610932f);
    for (int i=0; i<20; i++) {
      int  i1=indexData[i][0];
      int  i2=indexData[i][1];
      int  i3=indexData[i][2];
      Vector v1=new Vector(vdata[i1][0], vdata[i1][1], vdata[i1][2]);
      Vector v2=new Vector(vdata[i2][0], vdata[i2][1], vdata[i2][2]);
      Vector v3=new Vector(vdata[i3][0], vdata[i3][1], vdata[i3][2]); 
      v1=mat0.Mult(mat0, v1);
      v2=mat0.Mult(mat0, v2);
      v3=mat0.Mult(mat0, v3);
      Vector c1=NoiseVector(v1, planet.mass);
      Vector c2=NoiseVector(v2, planet.mass);
      Vector c3=NoiseVector(v3, planet.mass);
      int fracLevel=depth;
      Vector center=v1.Add(v1, v2).AsVector();
      center.Scale(0.5);
      center=center.Add(center, v2).AsVector();
      center.Scale(0.5);
      /*for (int j=0; j<fracLevel; j++) {
       c1=c1.Add(c1, new Vector(arcRandomAverage(planet.wind), arcRandomAverage(planet.land), arcRandomAverage(planet.rain))).AsVector();
       c2=c2.Add(c2, new Vector(arcRandomAverage(planet.wind), arcRandomAverage(planet.land), arcRandomAverage(planet.rain))).AsVector();
       c3=c3.Add(c3, new Vector(arcRandomAverage(planet.wind), arcRandomAverage(planet.land), arcRandomAverage(planet.rain))).AsVector();
       c1=c1.Add(c1, NoiseVector(v1, planet.mass)).AsVector();
       c2=c2.Add(c2, NoiseVector(v2, planet.mass)).AsVector();
       c3=c3.Add(c3, NoiseVector(v3, planet.mass)).AsVector();
       
       int rndLevel=depth+arcRandomInt(6);
       for (int k=0; k<rndLevel; k++) {
       Vector r1=c1.Add(c2, c3).AsVector();
       Vector r2=c2.Add(c1, c3).AsVector();
       Vector r3=c3.Add(c2, c1).AsVector();
       double t1=noise((float)r1.y, (float)r1.z);
       double t2=noise((float)r2.y, (float)r2.z);
       double t3=noise((float)r3.y, (float)r3.z);
       r1.Scale((AverageNoise(c1.x, c1.y, c1.z, planet.size)+T+Math.abs(v1.y*t1))/(c1.Length()));
       r2.Scale((AverageNoise(c2.x, c2.y, c2.z, planet.size)+T+Math.abs(v2.y*t2))/(c2.Length()));
       r3.Scale((AverageNoise(c3.x, c3.y, c3.z, planet.size)+T+Math.abs(v3.y*t3))/(c3.Length()));
       c1=c1.Add(c1, r1).AsVector();
       c2=c2.Add(c2, r2).AsVector();
       c3=c3.Add(c3, r3).AsVector();
       }
       }
       c1.Scale(0.25/fracLevel);
       c2.Scale(0.25/fracLevel);
       c3.Scale(0.25/fracLevel);*/
      Clamp(c1);
      Clamp(c2);
      Clamp(c3);
      //
      double h1=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*planet.energy*(1+Math.abs(v1.y));
      double h2=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*planet.energy*(1+Math.abs(v2.y));
      double h3=(arcRandomAverage(planet.rivers)-arcRandomAverage(planet.rain))*planet.energy*(1+Math.abs(v3.y));
      Vertex vx1=new Vertex(v1, v1.AsNormal(), c1, h1);
      Vertex vx2=new Vertex(v2, v2.AsNormal(), c2, h2);
      Vertex vx3=new Vertex(v3, v3.AsNormal(), c3, h3);
      vn[i][0]=vx1;
      vn[i][1]=vx2;
      vn[i][2]=vx3;
      vn[i][0].index=i1;
      vn[i][1].index=i2;
      vn[i][2].index=i3;
    }
    Vertex verts[]=new Vertex[12];
    for (int i=0; i<20; i++) {
      for (int j=0; j<3; j++) {
        if (verts[vn[i][j].index]==null) {
          verts[vn[i][j].index]=vn[i][j];
        }
      }
      triangle[i]=new PlanetTessa(verts[vn[i][0].index], verts[vn[i][1].index], verts[vn[i][2].index]);
    }
    ge=new GraphicsOrder[20]; 
    //list of all verts
    for (int i=0; i<20; i++) {
      triangle[i].SubDivide(triangle[i].vert[0], triangle[i].vert[1], triangle[i].vert[2], planet, depth);
      ge[i]=new GraphicsOrder();
      //triangle[i].GetLevel(ge[i].vertecies, 6);
    }
  }
  void Path(int loops) {
    for (int l=0; l<loops; l++) {
      int d=arcRandomInt(20);
      Tessa start=null;
      while (start==null) {
        start=ge[d].vertecies.get(arcRandomInt(ge[d].vertecies.size()));
        d=arcRandomInt(20);
      }
      Tessa next=start;
      Vertex last=null;
      ArrayList<Vertex>vis=new ArrayList<Vertex>(); 
      for (int i=0; i<loops; i++) {

        //println("Path "+i);
        int v=arcRandomInt(3);
        int a=arcRandomInt(next.vert[v].adjacent);
        Vertex adj=next.vert[v].adj[a];

        if (!vis.contains(adj) && adj!=null) {
          mole1.SingleCalc();
          adj.h+=mole1.Cos()*planet.land;
          last=adj;
          vis.add(last);
          if (adj.parent!=null)
            next=adj.parent;
          //int b=arcRandomInt(4);
          //next=next.tess[b];
        }
      }
      d=arcRandomInt(20);
      start=null;
      while (start==null) {
        start=ge[d].vertecies.get(arcRandomInt(ge[d].vertecies.size())); 
        d=arcRandomInt(20);
      }
      next=start;
      last=null;
      for (int i=0; i<loops; i++) {
        //println("Path "+i);
        int v=arcRandomInt(3);
        int a=arcRandomInt(next.vert[v].adjacent);
        Vertex adj=next.vert[v].adj[a];
        if (!vis.contains(adj) && adj!=null) {
          mole1.SingleCalc();
          adj.h-=mole1.Sin()*planet.water;
          last=adj;
          vis.add(last);
          if (adj.parent!=null)
            next=adj.parent;
          //int b=arcRandomInt(4);
          //next=next.tess[b];
        }
      }
    }
  }
  void Fracture(int loops) {
    for (int l=0; l<loops; l++) {
      int d=arcRandomInt(20);
      Tessa start=null;
      while (start==null) {
        start=ge[d].vertecies.get(arcRandomInt(ge[d].vertecies.size()));
        d=arcRandomInt(20);
      }
      Tessa next=start;
      Vertex last=null;
      ArrayList<Vertex>vis=new ArrayList<Vertex>(); 
      for (int i=0; i<loops; i++) {

        //println("Path "+i);
        int v=arcRandomInt(3);
        int a=arcRandomInt(next.vert[v].adjacent);
        Vertex adj=next.vert[v].adj[a];

        if (!vis.contains(adj) && adj!=null) {
          adj.h=(adj.h*adj.h)*0.5-adj.h*planet.land;
          last=adj;
          vis.add(last);
          if (adj.parent!=null)
            next=adj.parent;
          //int b=arcRandomInt(4);
          //next=next.tess[b];
        }
      }
      d=arcRandomInt(20);
      start=null;
      while (start==null) {
        start=ge[d].vertecies.get(arcRandomInt(ge[d].vertecies.size())); 
        d=arcRandomInt(20);
      }
      next=start;
      last=null;
      for (int i=0; i<loops; i++) {
        //println("Path "+i);
        int v=arcRandomInt(3);
        int a=arcRandomInt(next.vert[v].adjacent);
        Vertex adj=next.vert[v].adj[a];
        if (!vis.contains(adj) && adj!=null) {

          adj.h=(adj.h*adj.h)*0.5-adj.h*planet.water;
          last=adj;
          vis.add(last);
          if (adj.parent!=null)
            next=adj.parent;
          //int b=arcRandomInt(4);
          //next=next.tess[b];
        }
      }
    }
  }
  void Diffusion(int loops, int probability) {
    GraphicsOrder go=new GraphicsOrder();
    ArrayList<Integer> index=new ArrayList<Integer>(); 
    GraphicsOrder fixed=new GraphicsOrder(); 
    for (int i=0; i<20; i++) {
      for (int j=0; j<ge[i].vertecies.size(); j++) {
        if (arcRandomInt(probability)==0) {
          fixed.vertecies.add(ge[i].vertecies.get(j));
        } else {
          go.vertecies.add(ge[i].vertecies.get(j));
        }
      }
    }
    for (int j=0; j<go.vertecies.size(); j++) {
      index.add(new Integer(fixed.vertecies.size()+arcRandomInt(go.vertecies.size())));
    }
    for (int l=0; l<loops; l++) {
      //
      int ls=l+1;
      for (int i=0; i<go.vertecies.size(); i++) {
        int ix=index.get(i);
        ix+=arcRandomInt(go.vertecies.size())-arcRandomInt(go.vertecies.size());
        index.set(i, new Integer(ix));
        if (ix>0 && ix<fixed.vertecies.size()) {
          Tessa t1=fixed.vertecies.get(ix);
          int n=arcRandomInt(3);
          Tessa t3=t1.vert[n].adj[arcRandomInt(t1.vert[n].adjacent)].parent;
          if (!fixed.vertecies.contains(t3)) {
            fixed.vertecies.add(t3);
          }
          Tessa t2=go.vertecies.get(i);
          go.vertecies.remove(i);
          index.remove(i);
          double n1=noise(i*ls, ix*ls)*planet.land-noise(ix*ls, i*ls)*planet.water;
          double n2=n1*n1;
          for (int j=0; j<3; j++) {            
            t1.vert[j].h+=n2+n1*j;
            t3.vert[j].h+=t2.vert[j].h;
            t3.vert[j].h*=0.5;
            t2.vert[j].h-=n2-n1*j;
            t2.vert[j].h*=0.5;
          }
          //
        }
      }
    }
  }
  void SetVerts() {

    for (int i=0; i<20; i++) {
      ge[i]=new GraphicsOrder();
      triangle[i].GetLevel(ge[i].vertecies, depth);
    }
    Equality();
    //
    for (int i=0; i<20; i++) {
      triangle[i].Guide(triangle[i]);
    }
    for (int i=0; i<20; i++) {
      ge[i]=new GraphicsOrder();
      triangle[i].GetLevel(ge[i].vertecies, depth);
    }
    Adjacancy();
    //
  }
  void Equality() {
    //merge duplicate verts 
    int eq=0;
    for (int i=0; i<20; i++) {
      println("Equity "+i+" Eq "+eq+" Size "+ge[i].vertecies.size());
      eq=0;
      //if(i!=15){//no equality
      for (int j=0; j<20; j++) {
        int ind=0;
        //if(i!=j){
        for (int ii=0; ii<3; ii++) {

          for (int jj=0; jj<3; jj++) {
            if (triangle[i].vert[ii].index==triangle[j].vert[jj].index) {
              ind++;
            }
          }
        }
        boolean testOk=ind>1;
        if (testOk) {
          for (int k1=0; k1<ge[i].vertecies.size (); k1++) {
            Tessa test=ge[i].vertecies.get(k1);
            for (int k2=0; k2<ge[j].vertecies.size (); k2++) {
              Tessa against=ge[j].vertecies.get(k2);
              if (test!=against) {
                eq+=against.Equity(test);
              }
              //triangle[j].Equality(triangle[i], 7);
            }
          }
        }
        //}
      }
    }
  }
  void Adjacancy() {
    int ja=0;
    for (int i=0; i<20; i++) {
      println("Adjacancy "+i+" adj "+ja);
      ja=0;
      for (int j=0; j<20; j++) {
        //if(i!=j){
        for (int ai=0; ai<ge[i].vertecies.size (); ai++) {
          /*double n6[][]= {
           {
           MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP
           }
           , 
           {
           MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP
           }
           , 
           {
           MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP, MAX_MAP
           }
           };*/
          Tessa test=ge[i].vertecies.get(ai);
          int ind=0;
          for (int ii=0; ii<3; ii++) {

            for (int jj=0; jj<3; jj++) {
              if (triangle[i].vert[ii].index==triangle[j].vert[jj].index) {
                ind++;
              }
            }
          }
          boolean testOk=ind>1;
          if (testOk) {
            double nc=MAX_MAP;
            //int n[]=new int[3];
            for (int aj=0; aj<ge[j].vertecies.size (); aj++) {

              Tessa adjTest=ge[j].vertecies.get(aj);
              //if(adjTest!=test){
              for (int v1=0; v1<3; v1++) {
                for (int v2=0; v2<3; v2++) {
                  if (test.vert[v1]==adjTest.vert[v2]) {
                    if (test.vert[v1].adjacent>0) {
                      int sofar=test.vert[v1].adjacent;
                      for (int v3=0; v3<3; v3++) {
                        boolean ok=true;
                        for (int na=0; na<sofar; na++) {
                          if (v3==v2 || test.vert[v1].adj[na]==adjTest.vert[v3]) {
                            ok=false;
                          }
                        }
                        if (ok) { 
                          test.vert[v1].adj[test.vert[v1].adjacent]=adjTest.vert[v3];
                          test.vert[v1].adjacent++;
                          ja++;
                        }
                      }
                    } else {
                      for (int v3=0; v3<3; v3++) {
                        if (v3!=v2 && test.vert[v1].adj[test.vert[v1].adjacent]!=adjTest.vert[v3]) {
                          test.vert[v1].adj[test.vert[v1].adjacent]=adjTest.vert[v3];
                          test.vert[v1].adjacent++;
                          ja++;
                        }
                      }
                    }
                  }
                }
                //}
              }
              /*   
               Cartesian cdis=test.center.Sub(test.center, adjTest.center);
               double cd=cdis.LengthSquared();
               if (cd<=nc) {
               nc=cd;
               for (int n=0; n<6; n++) {
               for (int v2=0; v2<3; v2++) {
               for (int v1=0; v1<3; v1++) {
               //
               if (test.vert[v1]!=adjTest.vert[v2]) {
               Cartesian dis=test.vert[v1].v.Sub(test.vert[v1].v, adjTest.vert[v2].v);
               double d=dis.LengthSquared();
               if (d<n6[v1][n]) {
               n6[v1][n]=d;
               adj[v1][n]=adjTest.vert[v2];
               
               }
               }
               }
               //
               }
               }
               */
            }
            //
          }
          //}
          //
        }
      }
    }
  }
  void Rivers(int w) {
    for (int n=0; n<w; n++) {
      int n8=n%10;
      for (int i=0; i<20; i++) {

        //if(triangle[i].vert[0].h>0 ||triangle[i].vert[1].h>0 ||triangle[i].vert[2].h>0){
        triangle[i].RiverLevel(n8, 0);
        //}
      }
    }
  }
  void LandFalls(int w) {
    double im=1.0/planet.mass;
    for (int n=0; n<w; n++) {
      int n8=n%10;
      for (int i=0; i<20; i++) {

        triangle[i].LandFalls(n8, im, depth);
      }
    }
  }
  void River() {

    for (int n=1; n<10; n++) {  
      for (int i=0; i<20; i++) {
        triangle[i].Riverer(planet, n);
      }
    }
  }
  void Merge(int m) { 
    for (int i=0; i<m; i++) {

      int i1=arcRandomInt(20);
      int i2=arcRandomInt(20);
      if (i1!=i2)triangle[i1].Merge(triangle[i2]);
    }
  }
  void Blur() {

    for (int n=1; n<10; n++) {  
      for (int i=0; i<20; i++) {

        for (int k=0; k<ge[i].vertecies.size (); k++) {
          Tessa ts=ge[i].vertecies.get(k);
          triangle[i].Blur(n, 0, ts);
        }
      }
    }
  }
  void Contour() {
    for (int n=1; n<10; n++) {
      for (int i=0; i<20; i++) {
        triangle[i].ContourLevel(planet, n);
      }
    }
  }
  void Transfer() {

    for (int n=1; n<10; n++) {
      for (int i=0; i<20; i++) {
        triangle[i].Transfer(n);
      }
    }
  }
  void Blend() {

    for (int i=0; i<20; i++) {
      triangle[i].Blender(planet, depth);
    }
  }
  void Flatten() {

    for (int i=0; i<20; i++) {
      triangle[i].Flattener(depth);
    }
  }
  void Colorize(Vector stain) {

    for (int i=0; i<20; i++) {
      triangle[i].ColorizeLevel(planet, stain, depth);
    }
  }
  void AddSphere() {

    for (int i=0; i<20; i++) {

      triangle[i].Add(triangleBuffer[i]);
    }
  }
  void ScaleSphere(double scale) {

    for (int i=0; i<20; i++) {

      triangle[i].Scale(triangle[i], scale);
    }
  }
  ArrayList<Tessa> GetVerts(Vector cam, int depth) {
    ArrayList<Tessa> list=new ArrayList<Tessa>();
    for (int i=0; i<20; i++) {

      triangle[i].GetLevel(list, depth);
    }
    return list;
  }
}
