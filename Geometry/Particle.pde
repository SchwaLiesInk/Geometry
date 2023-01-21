static int b=150;
static int PARTICLES=100;
class ParticleDisplacer extends Point implements VectorStateFunction {
  ParticleDisplacer() {
    r=10;
  }
    Vector fv(double x, double y, double z, double w) {
    return new Vector(x*r*Math.sin(w), y, z*r*Math.cos(w));
  }
}
class Particle implements OneStateFunction {
  double at=0;
  double from=0;
  double to=0;
  double point=0;
  //
  Point local;
  Vector direct;
  Vector velocity;
  //
  double startingRadius=5;
  double splash=3;
  double angle=0;
  double decay=0.1;
  double gravity=9.8;
  double spread=1;
  double travelled=0;
  double travels=b;
  double energy=3;
  double mass=3;
  double scale=100;
  //
  double a=1;//potential energy start (symetrical)
  double normal =-10;//normalizing constant(fire) chaos if negative -e(air preasure)
  //
  ParticleIterator iterate=new ParticleIterator();
  ParticleResult result;
  ParticleDisplacer placer=new ParticleDisplacer();
  NumberSequance ns;
  IntegralOperation in;
  //
  Particle() {
    local=new Point(0, 0, 0, startingRadius);
    direct=new Vector(arcRandom()-arcRandom(), arcRandom()-arcRandom(), arcRandom()-arcRandom());
    direct.Normalize();
    velocity=new Vector();
    ns= new NumberSequance(new Double[] {
      (double)0.0
      }
      );
    result=new ParticleResult(1+(1+arcRandomInt(4))*arcRandom(), 1+(1+arcRandomInt(5))*arcRandom());

    in = new IntegralOperation(a, b+0.2, normal);
    in.SetVariables(ns);
    in.OperateUsing(this);//scale function
    in.iterate=iterate;//scaler
    in.resultant=result;
    splash=1;
  }
  double Move(double t, double lim, double dampen, double spin) {
    Matrix rot=new Matrix();
    rot=rot.RotateY(spin*t);
    angle+=spin*t;
    direct=rot.Mult(rot, direct);
    //
    double acc=ac(t);
    double acv=scale*acc/mass*t;
    Vector i=direct.Scaled(acv);
    Vector g=new Vector(0, scale, 0).Scaled(gravity*t);
    Vector p=new Vector(0, -scale, 0).Scaled(acv*splash*t);
    velocity=velocity.Add(velocity, i).AsVector();
    velocity=velocity.Add(velocity, g).AsVector();
    velocity=velocity.Add(velocity, p).AsVector();
    if (placer!=null) {
      velocity=velocity.Add(velocity, placer.fv(acv, i.y, acv, angle)).AsVector();
    }
    velocity=velocity.Scaled(1.0-dampen);
    Vector displace=velocity.Scaled(t);
    //
    local=(local.Add(local, displace)).AsPoint(local.r-local.r*decay*t);
    double r=local.Length();
    double v=velocity.Length();
    splash+=spread/(gravity*mass/local.r);//
    //
    travelled+=v*t;
    if (travelled>travels || r>lim || local.r<1.0) {
      direct=new Vector(arcRandom()-arcRandom(), arcRandom()-arcRandom(), arcRandom()-arcRandom());
      direct.Normalize();
      velocity=new Vector();
      local.x=0;
      local.y=0;
      local.z=0;
      local.r=startingRadius;
      in.ResetState();
      splash=1;
      travelled=0;
    }
    return r;
  }
    double f(double x) {
    return (1.0/x-x*(1.0-1.0/b));
  }
  double ac(double t) {

    if (point==0) {
      in.LoopSum(t);
      to=in.results.First();
      at=from;
    }
    if (point<1.0+t) {
      at=lerp(point, from, to);
      point+=t;
    } else {
      from=to;
      point=0;
    }
    if (in.results.First()!=null) {
      return in.results.First()*energy;
    } else {
      return 0;
    }
  }
  void Print() {

    int i=0;
    do {
      double t=0.1;//(double)(i+1)/(b+(1+i));
      in.LoopSum(t);//energy system
      double sys=(double)(i/b);
      System.out.println(sys+" Result "+in.a+" << \t" + in.results.ToString(10));
      if (point==0) {
        to=in.results.First();
        at=from;
      }
      while (point<1.0+t) {
        //System.out.println(" Operation "+point+" Lerp "+at);
        at=lerp(point, from, to);
        point+=t;
      }
      from=to;
      point=0;
      /*if(in.results.First()<0.1){
       in.a=Random.arcRandom();//add fuel
       in.c=-Random.arcRandom();
       }
       //in.c-=t;//fire entropy
       //System.out.println(sys+" Rate Of Change "+in.c+" >> " + in.changes.ToString(10));
       */
      i++;
    } while (!in.Finalize ());
  }
}
class ParticleIterator implements OneStateFunction {

    double f(double x) {
    return (x/b);
  }
}
class ParticleResult implements OneStateFunction {
  double dim=4.0;
  double entropy=4.5;
  ParticleResult(double d, double en) {
    dim=d;
    entropy=en;
  }
    double f(double x) {
    return entropy/(x*dim/3.0);
  }
}
