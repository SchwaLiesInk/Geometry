class Molecula {
  double cosine;
  double sine;
  double startAngle;
  double currentAngle;
  double spin;
  ///
  Molecula() {
    startAngle=Math.PI*arcRandom()/arcRandomInt(3600);
    currentAngle=PI2*arcRandom();
    spin=startAngle*arcRandom();
    cosine=Math.cos(spin);
    sine=Math.sin(spin);
  }
  Molecula(double start, double current, double spin) {
    startAngle=start;
    currentAngle=current;
    this.spin=startAngle*spin;
    cosine=Math.cos(spin);
    sine=Math.sin(spin);
  }
  String ToString() {
    return "Molecula angles\tstart"+startAngle+"\tcurrent"+currentAngle+"\tspin"+spin+"\ttrig cos"+cosine+"sine"+cosine;
  }
  double Angle()
  {
    return this.currentAngle;
  }
  double Spin()
  {
    return this.spin;
  }
  double Cos()
  {
    return cosine;
  }
  double Sin()
  {
    return sine;
  }
  void SingleCalc() {
    spin+=currentAngle;
    cosine=Math.cos(spin);
    sine=Math.sin(spin);
    currentAngle+=startAngle;
  }
  void SingleCalc(double scale) {
    spin+=currentAngle;
    cosine=Math.cos(spin);
    sine=Math.sin(spin);
    currentAngle+=startAngle*scale;
  }
  void Calc() {
    spin+=currentAngle;
    cosine=Math.cos(spin);
    spin+=currentAngle;
    sine=Math.sin(spin);
    currentAngle+=startAngle;
  }
  void Calc(double scale) {
    spin+=currentAngle;
    cosine=Math.cos(spin);
    spin+=currentAngle;
    sine=Math.sin(spin);
    currentAngle+=startAngle*scale;
  }
  void LineCalc() {
    spin+=currentAngle;
    spin/=startAngle;
    cosine=Math.cos(spin);
    spin+=currentAngle;
    spin/=startAngle;
    sine=Math.sin(spin);
    currentAngle+=startAngle;
  }
  void LineCalc(double scale) {
    spin+=currentAngle;
    spin/=startAngle;
    cosine=Math.cos(spin);
    spin+=currentAngle;
    spin/=startAngle;
    sine=Math.sin(spin);
    currentAngle+=startAngle*scale;
  }
}
double lerp(double  t, double v1, double v2) { 
  return (1.0f - t) * v1 + t * v2;
}

Random aRandom=new Random();
Molecula mole1;
Molecula mole2;
Molecula mole3;
//Particle aParticle;
//float particleLoop=1.0f;
void arcRandomSeed(long s) {
  aRandom=new Random(s);
  mole1=new Molecula();
  mole2=new Molecula();
  mole3=new Molecula();
  //aParticle=new Particle();
}
float arcRandom(float d) {
  //particleLoop=(float)aParticle.Move(0.005,1.0,0.01,100);
  // println("L:"+particleLoop+"XZ:"+Math.abs(aParticle.local.x+aParticle.local.z)*5);
  return aRandom.nextFloat()*d;
}
float arcRandom() {
  //particleLoop=(float)aParticle.Move(0.005, 1.0, 0.01, 100);
  //println("L:"+particleLoop+"XZ:"+Math.abs(aParticle.local.x+aParticle.local.z)*5);
  //float f=(float)Math.abs(aParticle.local.x+aParticle.local.z)*5;
  return aRandom.nextFloat();
}
int arcRandomInt(int p) {
  //particleLoop=(float)aParticle.Move(0.005,1.0,0.01,100);
  //println("L:"+particleLoop+"XZ:"+Math.abs(aParticle.local.x+aParticle.local.z)*5);
  //float f=(float)Math.abs(aParticle.local.x+aParticle.local.z)*5;
  return aRandom.nextInt(p);
}
boolean arcRandomBoolean() {
  //particleLoop=(float)aParticle.Move(0.005,1.0,0.01,100);
  //println("L:"+particleLoop+"XZ:"+Math.abs(aParticle.local.x+aParticle.local.z)*5);
  //float f=(float)Math.abs(aParticle.local.x+aParticle.local.z)*5;
  return aRandom.nextBoolean();
}
int arcRandomInt(int n, int to) {
  return n + arcRandomInt(1 + to - n);
}
double T=0.23456789123456789123456789;
float S=2.0f;
double arcRandomAverage(double rsc) {
  mole1.SingleCalc();
  S=(float)(mole1.Cos()*mole1.Sin());
  //return noise((1+S)+(float)mole1.Cos()*(1+S), (1+S)+(float)mole1.Sin()*(1+S))*arcRandom()*S*scale;
  return (arcRandom()-arcRandom())*arcRandom()*S*rsc;
}
Vector NoiseVector(Vector v1, double scale) {
  mole2.SingleCalc();
  mole1.SingleCalc();
  S=(float)(arcRandomAverage((mole1.Cos()+mole1.Sin()))*scale);
  return new Vector(Noise((float)((1+v1.x)*S), (float)(1+mole2.Cos())*S)*S, Noise((float)(1+v1.y)*S, (float)(1+mole2.Sin())*S)*S, Noise((float)(1+v1.z)*S, (float)((2+mole2.Cos()+mole2.Sin())*S))*S);
}
double AverageNoise(double x, double y, double z, double scale) {
  mole3.SingleCalc();
  mole1.SingleCalc();
  S=(float)arcRandomAverage((mole1.Cos()+mole1.Sin()));
  return Math.abs(Noise((float)((1+x)*S), (float)((1+mole3.Cos())*S))-Noise((float)((1+y)*S), (float)((1+mole3.Sin())*S)))*Noise((float)((1+z)*S), (float)(((2+mole3.Cos()+mole3.Sin())*S)))*scale*S;
}

double Noise(double x, double y) {
  double ir=0.456789123/Math.sqrt(1+x*x+y*y);
  return (arcRandom()*(0.1+Math.abs(x))*ir+arcRandom()*(0.1+Math.abs(y))*ir);
}
double Noise(double x) {
  double ir=0.987654321/Math.sqrt(1+x*x);
  return arcRandom()*(0.1+Math.abs(x))*ir;
}

int dice(int prob) {
  return 1 + arcRandomInt(1 + prob);
}

int dice(int n, int prob) {
  int dn = 0;
  for (int i = 0; i < n; i++) {
    dn += 1 + arcRandomInt(1 + prob);
  }
  return dn;
}


interface OneStateFunction {
  public abstract double f(double x);
}
interface TwoStateFunction {//reducer
  public abstract double f2(double x, double y);
}
interface ThreeStateFunction {
  public abstract double f3(double x, double y, double z);
}
interface VectorStateFunction {
  public abstract Vector fv(double x, double y, double z, double w);
}
interface FourStateFunction {
  public abstract double f4(Vector v);
}
interface State3DFunction {
  public abstract Vector f3D(Matrix m);
}
interface MultiStateFunction {
  public abstract double g(RingList<Double>xyz);
}
interface OneSumFunction {
  double h(OneStateFunction fn);
}
interface TwoSumFunction {
  double h2(TwoStateFunction fn);
}
interface ThreeSumFunction {
  double h3(ThreeStateFunction fn);
}
interface FourSumFunction {
  Vector h4(FourStateFunction fn);
}
interface MultiSumFunction {
  RingList<Double> i(MultiStateFunction fn);
}

interface Sum3DFunction {
  public abstract Matrix f3D(State3DFunction m);
}
class Maths {

  /*public static NumberSequance TestFiber() {
   
   NumberSequance ns = new NumberSequance(NumberSequance.OneTo9, (nb) -> {
   return nb + nb + 1;
   });
   System.out.println(ns.ToString(10));
   ns.Calculate();
   System.out.println(ns.ToString(10));
   return ns;
   }
   
   public static NumberSequance TestFiber(int n) {
   
   NumberSequance ns = new NumberSequance();
   ns.Generate((OneStateFunction) (nb) -> nb + nb + 1, 1, n, 1);
   System.out.println(ns.ToString(10));
   return ns;
   }
   
   public static NumberSequance TestSine(int n, double seed) {
   
   NumberSequance ns = new NumberSequance();
   ns.GenerateLoop(n, seed, (OneStateFunction) (nb) -> Math.sin(nb) * Linear.PI2);
   System.out.println(ns.ToString(10));
   return ns;
   }
   
   public static void TestTrig(NumberSequance ns, OneStateFunction trig) {
   
   ns = ns.Divide(new NumberSequance(NumberSequance.TwoTo10), (n) -> n * Math.PI);
   System.out.println(ns.ToString(10));
   ns = ns.Operation((OneStateFunction) trig);
   
   System.out.println(ns.ToString(10));
   }
   
   
   
   public static NumberSequance TestShuffle() {
   
   NumberSequance ns = new NumberSequance(NumberSequance.TwoTo10);
   System.out.println(ns.ToString(10));
   ns = ns.Shuffle(0);
   System.out.println(ns.ToString(10));
   return ns;
   }*/
}
final static double OneTo9[] = new double[]{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
final static double OneTo10PI[] = new double[]{Math.PI, 2.0 * Math.PI, 3.0 * Math.PI, 4.0 * Math.PI, 5.0 * Math.PI, 6.0 * Math.PI, 7.0 * Math.PI, 8.0 * Math.PI, 9.0 * Math.PI, 10.0 * Math.PI};
final static double TwoTo10[] = new double[]{2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

class NumberSequance extends RingList<Double> implements OneStateFunction {

  int n;
  int a;
  int b;
  int step;
  OneStateFunction function;
  Double currentValue;
  Node<Double> currentNode;

  NumberSequance() {
    super();
  }

  NumberSequance(double v) {
    super();
    Append(v);
  }

  NumberSequance(double v1, double v2) {
    super();
    Append(v1);
    Append(v2);
  }

  NumberSequance(double v1[], double v2[]) {

    super();
    for (int i = 0; i < v1.length && i < v2.length; i++) {
      Append(v1[i]);
      Append(v2[i]);
    }
  }

  NumberSequance(double v1[], double v2[], double v3[]) {

    super();
    for (int i = 0; i < v1.length && i < v2.length && i < v3.length; i++) {
      Append(v1[i]);
      Append(v2[i]);
      Append(v3[i]);
    }
  }

  NumberSequance(double v1[], double v2[], double v3[], double v4[]) {

    super();
    for (int i = 0; i < v1.length && i < v2.length && i < v3.length && i < v4.length; i++) {
      Append(v1[i]);
      Append(v2[i]);
      Append(v3[i]);
      Append(v4[i]);
    }
  }

  NumberSequance(double v1, double v2, double v3) {
    super();
    Append(v1);
    Append(v2);
    Append(v3);
  }

  NumberSequance(Point p) {
    super();
    Append(p.x);
    Append(p.y);
    Append(p.z);
    Append(p.r);
  }

  NumberSequance(Point p[]) {
    super();
    for (int i = 0; i < p.length; i++) {
      Append(p[i].x);
      Append(p[i].y);
      Append(p[i].z);
      Append(p[i].r);
    }
  }

  NumberSequance(Vector v) {
    super();
    Append(v.x);
    Append(v.y);
    Append(v.z);
    Append(v.w);
  }

  NumberSequance(Vector v[]) {
    super();
    for (int i = 0; i < v.length; i++) {
      Append(v[i].x);
      Append(v[i].y);
      Append(v[i].z);
      Append(v[i].w);
    }
  }

  NumberSequance(Cartesian v[], OneStateFunction com) {
    super();
    for (int i = 0; i < v.length; i++) {
      Append(com.f(v[i].x));
      Append(com.f(v[i].y));
      Append(com.f(v[i].z));
    }
  }

  NumberSequance(Cartesian v1[], Cartesian v2[], TwoStateFunction com) {
    super();
    for (int i = 0; i < v1.length && i < v2.length; i++) {
      Append(com.f2(v1[i].x, v2[i].x));
      Append(com.f2(v1[i].y, v2[i].y));
      Append(com.f2(v1[i].z, v2[i].z));
    }
  }

  NumberSequance(Cartesian v1[], Cartesian v2[], Cartesian v3[], ThreeStateFunction com) {
    super();
    for (int i = 0; i < v1.length && i < v2.length && i < v3.length; i++) {
      Append(com.f3(v1[i].x, v2[i].x, v3[i].x));
      Append(com.f3(v1[i].y, v2[i].y, v3[i].y));
      Append(com.f3(v1[i].z, v2[i].z, v3[i].z));
    }
  }

  NumberSequance(Cartesian v[], OneStateFunction comx, OneStateFunction comy, OneStateFunction comz) {
    super();
    for (int i = 0; i < v.length; i++) {
      Append(comx.f(v[i].x));
      Append(comy.f(v[i].y));
      Append(comz.f(v[i].z));
    }
  }

  NumberSequance(Cartesian v1[], Cartesian v2[], TwoStateFunction comx, TwoStateFunction comy, TwoStateFunction comz) {
    super();
    for (int i = 0; i < v1.length && i < v2.length; i++) {
      Append(comx.f2(v1[i].x, v2[i].x));
      Append(comy.f2(v1[i].y, v2[i].y));
      Append(comz.f2(v1[i].z, v2[i].z));
    }
  }

  NumberSequance(Cartesian v1[], Cartesian v2[], Cartesian v3[], ThreeStateFunction comx, ThreeStateFunction comy, ThreeStateFunction comz) {
    super();
    for (int i = 0; i < v1.length && i < v2.length && i < v3.length; i++) {
      Append(comx.f3(v1[i].x, v2[i].x, v3[i].x));
      Append(comy.f3(v1[i].y, v2[i].y, v3[i].y));
      Append(comz.f3(v1[i].z, v2[i].z, v3[i].z));
    }
  }

  NumberSequance(Cartesian v) {
    super();
    Append(v.x);
    Append(v.y);
    Append(v.z);
  }

  NumberSequance(Cartesian v[]) {
    super();
    for (int i = 0; i < v.length; i++) {
      Append(v[i].x);
      Append(v[i].y);
      Append(v[i].z);
    }
  }

  NumberSequance(Double v[]) {
    super(v);
  }

  NumberSequance(RingList<Double> v) {
    super(v);
  }

  NumberSequance(Collection<Double> v) {
    super(v);
  }

  NumberSequance(double v[]) {
    super();
    for (int i = 0; i < v.length; i++) {
      Append(v[i]);
    }
  }

  NumberSequance(OneStateFunction fn) {
    super();
    this.function = fn;
  }

  NumberSequance(Double v[], OneStateFunction fn) {
    super(v);
    this.function = fn;
  }

  NumberSequance(RingList<Double> v, OneStateFunction fn) {
    super(v);
    this.function = fn;
  }

  NumberSequance(Collection<Double> v, OneStateFunction fn) {
    super(v);
    this.function = fn;
  }

  NumberSequance(double v[], OneStateFunction fn) {
    super();
    for (int i = 0; i < v.length; i++) {
      Append(v[i]);
    }
    this.function = fn;
  }

  NumberSequance Shuffle(int times) {
    times = Math.abs(times);
    NumberSequance odd = this.Odd();
    NumberSequance even = this.Even();
    //System.out.println("odd "+odd.length);
    //System.out.println("even "+even.length);
    NumberSequance result = new NumberSequance();
    result = even.Combined(odd);
    int i = 0;

    for (int t = 0; t < times; t++) {
      i = arcRandomInt(24);
      switch (i) {
      case 0: 
        {
          even.Reverse();
          result = even.Combined(odd);

          break;
        }
      case 1: 
        {
          odd.Reverse();
          result = odd.Combined(even);

          break;
        }
      case 2: 
        {
          odd.Reverse();
          result = even.Combined(odd);

          break;
        }
      case 3: 
        {
          even.Reverse();
          result = odd.Combined(even);

          break;
        }
      case 4: 
        {
          result = odd.Combined(even).Reversed();
          break;
        }
      case 5: 
        {
          result = even.Combined(odd).Reversed();
          break;
        }
      case 6: 
        {
          even.Reverse();
          result = odd.Inserted(even, arcRandomInt(odd.length));//cut

          break;
        }
      case 7: 
        {
          odd.Reverse();
          result = even.Inserted(odd, arcRandomInt(even.length));//cut

          break;
        }
      case 8: 
        {
          result = odd.Inserted(even, arcRandomInt(odd.length));//cut

          break;
        }
      case 9: 
        {
          result = even.Inserted(odd, arcRandomInt(even.length));//cut

          break;
        }
      case 10: 
        {
          result = odd.Inserted(even, arcRandomInt(odd.length)).Reversed();//cut

          break;
        }
      case 11: 
        {
          result = even.Inserted(odd, arcRandomInt(even.length)).Reversed();//cut

          break;
        }

      case 12: 
        {
          result = even.Intergrated(odd).Reversed();//cut
          break;
        }
      case 13: 
        {
          result = odd.Intergrated(even).Reversed();//cut
          break;
        }
      case 14: 
        {
          result = even.Intergrated(odd.Reversed());//cut
          break;
        }
      case 15: 
        {
          result = odd.Intergrated(even.Reversed());//cut
          break;
        }
      case 16: 
        {
          result = even.Intergrated(odd);//cut
          break;
        }
      case 17: 
        {
          result = odd.Intergrated(even);//cut
          break;
        }

      case 18: 
        {
          even.Reverse();
          result = odd.Intergrated(even, arcRandomInt(odd.length));//cut

          break;
        }
      case 19: 
        {
          odd.Reverse();
          result = even.Intergrated(odd, arcRandomInt(even.length));//cut

          break;
        }
      case 20: 
        {
          result = odd.Intergrated(even, arcRandomInt(odd.length));//cut

          break;
        }
      case 21: 
        {
          result = even.Intergrated(odd, arcRandomInt(even.length));//cut

          break;
        }
      case 22: 
        {
          result = odd.Intergrated(even, arcRandomInt(odd.length)).Reversed();//cut

          break;
        }
      case 23: 
        {
          result = even.Intergrated(odd, arcRandomInt(even.length)).Reversed();//cut

          break;
        }
      }
      i = arcRandomInt(2);
      if (i == 0) {
        odd = result.Even();
        even = result.Odd();
      } else {
        odd = result.Odd();
        even = result.Even();
      }
    }
    i = arcRandomInt(18);
    switch (i) {
    case 0: 
      {
        result = even.Inserted(odd, arcRandomInt(even.length)).Reversed();//cut
        break;
      }
    case 1: 
      {
        result = odd.Inserted(even, arcRandomInt(odd.length)).Reversed();//cut
        break;
      }
    case 2: 
      {
        result = even.Inserted(odd.Reversed(), arcRandomInt(even.length));//cut
        break;
      }
    case 3: 
      {
        result = odd.Inserted(even.Reversed(), arcRandomInt(odd.length));//cut
        break;
      }
    case 4: 
      {
        result = even.Inserted(odd, arcRandomInt(even.length));//cut
        break;
      }
    case 5: 
      {
        result = odd.Inserted(even, arcRandomInt(odd.length));//cut
        break;
      }
    case 6: 
      {
        result = even.Intergrated(odd).Reversed();//cut
        break;
      }
    case 7: 
      {
        result = odd.Intergrated(even).Reversed();//cut
        break;
      }
    case 8: 
      {
        result = even.Intergrated(odd.Reversed());//cut
        break;
      }
    case 9: 
      {
        result = odd.Intergrated(even.Reversed());//cut
        break;
      }
    case 10: 
      {
        result = even.Intergrated(odd);//cut
        break;
      }
    case 11: 
      {
        result = odd.Intergrated(even);//cut
        break;
      }

    case 12: 
      {
        even.Reverse();
        result = odd.Intergrated(even, arcRandomInt(odd.length));//cut

        break;
      }
    case 13: 
      {
        odd.Reverse();
        result = even.Intergrated(odd, arcRandomInt(even.length));//cut

        break;
      }
    case 14: 
      {
        result = odd.Intergrated(even, arcRandomInt(odd.length));//cut

        break;
      }
    case 15: 
      {
        result = even.Intergrated(odd, arcRandomInt(even.length));//cut

        break;
      }
    case 16: 
      {
        result = odd.Intergrated(even, arcRandomInt(odd.length)).Reversed();//cut

        break;
      }
    case 17: 
      {
        result = even.Intergrated(odd, arcRandomInt(even.length)).Reversed();//cut

        break;
      }
    }
    i = arcRandomInt(2);
    if (i == 0) {
      odd = result.Even();
      even = result.Odd();
    } else {
      odd = result.Odd();
      even = result.Even();
    }
    return even.Inserted(odd, arcRandomInt(even.length));//final cut
  }

  NumberSequance Reverse() {
    Stack<Double> st = new Stack();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      st.Push(d.data);
    }
    this.Clear();
    this.PasteAll(st.list);
    return this;
  }

  NumberSequance Reversed() {
    Stack<Double> st = new Stack();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      st.Push(d.data);
    }
    NumberSequance r = new NumberSequance();
    r.PasteAll(st.list);
    return r;
  }

  NumberSequance Clone() {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      r.Append(d.data);
    }
    return r;
  }

  NumberSequance Odd() {
    //zero counted as odd
    NumberSequance r = new NumberSequance();
    Node<Double> d = this.Start();
    for (; d.data != null && d.next.data != null; d = d.next.next) {
      r.Append(d.data);
    }
    if (d.data != null) {
      r.Append(d.data);
    }
    return r;
  }

  NumberSequance Even() {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null; d = d.next.next) {
      r.Append(d.next.data);
    }
    return r;
  }

  NumberSequance FunctionStep(OneStateFunction fn, int start, int step, boolean fill) {
    NumberSequance r = new NumberSequance();
    int i = 0;
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      if (start <= 0) {
        i++;
        if (i == step) {
          r.Append(fn.f(d.data));
          i = 0;
        } else if (fill) {
          r.Append(d.data);
        }
      } else {
        start--;
      }
    }
    return r;
  }

  //same size
  NumberSequance Function2Step(TwoStateFunction fn, int start, int step, boolean fill) {
    NumberSequance r = new NumberSequance();
    int i = 0;
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null; d = d.next) {
      if (start <= 0) {
        i++;
        if (i == step) {
          r.Append(fn.f2(d.data, d.next.data));
          i = 0;
        } else if (fill) {
          r.Append(d.data);
        }
      } else {
        start--;
      }
    }
    return r;
  }

  //same size
  NumberSequance Function3Step(ThreeStateFunction fn, int start, int step, boolean fill) {
    NumberSequance r = new NumberSequance();
    int i = 0;
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null && d.next.next.data != null; d = d.next) {
      if (start <= 0) {
        i++;
        if (i == step) {
          r.Append(fn.f3(d.data, d.next.data, d.next.next.data));
          i = 0;
        } else if (fill) {
          r.Append(d.data);
        }
      } else {
        start--;
      }
    }
    return r;
  }

  //same size
  NumberSequance Function(OneStateFunction fn) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      r.Append(fn.f(d.data));
    }
    return r;
  }

  //reduces
  NumberSequance Function(TwoStateFunction fn) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null; d = d.next.next) {
      r.Append(fn.f2(d.data, d.next.data));
    }
    return r;
  }

  //reduces
  NumberSequance Function(ThreeStateFunction fn) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null && d.next.next.data != null; d = d.next.next.next) {
      r.Append(fn.f3(d.data, d.next.data, d.next.next.data));
    }
    return r;
  }

  //expands
  NumberSequance Functions(OneStateFunction fn[]) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      for (int i = 0; i < fn.length; i++) {
        r.Append(fn[i].f(d.data));
      }
    }
    return r;
  }

  //same size/expand
  NumberSequance Functions(TwoStateFunction fn[]) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null; d = d.next.next) {
      for (int i = 0; i < fn.length; i++) {
        r.Append(fn[i].f2(d.data, d.next.data));
      }
    }
    return r;
  }

  //reduces/expands
  NumberSequance Functions(ThreeStateFunction fn[]) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null && d.next.next.data != null; d = d.next.next.next) {
      for (int i = 0; i < fn.length; i++) {
        r.Append(fn[i].f3(d.data, d.next.data, d.next.next.data));
      }
    }
    return r;
  }

  //x,y to y,x
  NumberSequance Flip2() {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null; d = d.next.next) {
      r.Append(d.next.data);
      r.Append(d.data);
    }
    return r;
  }

  //x,y,z to z,y,x
  NumberSequance Flip3() {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null && d.next.data != null && d.next.next.data != null; d = d.next.next.next) {
      r.Append(d.next.next.data);
      r.Append(d.next.data);
      r.Append(d.data);
    }
    return r;
  }

  //places at end
  NumberSequance Join(NumberSequance seq) {
    NumberSequance r = new NumberSequance();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      r.Append(new Double(d.data));
    }
    for (Node<Double> d = seq.Start(); d.data != null; d = d.next) {
      r.Append(new Double(d.data));
    }
    return r;
  }

  //interweaves
  NumberSequance Combined(NumberSequance seq) {
    NumberSequance r = new NumberSequance();
    Node<Double> d1 = this.Start();
    Node<Double> d2 = seq.Start();
    for (; d1.data != null && d2.data != null; d1 = d1.next, d2 = d2.next) {
      r.Append(new Double(d1.data));
      r.Append(new Double(d2.data));
    }
    if (d1.data != null) {
      for (; d1.data != null; d1 = d1.next) {
        r.Append(new Double(d1.data));
      }
    } else if (d2.data != null) {
      for (; d2.data != null; d2 = d2.next) {
        r.Append(new Double(d2.data));
      }
    }
    return r;
  }

  //interweaves
  NumberSequance Combined(NumberSequance seq, int step) {
    NumberSequance r = new NumberSequance();
    Node<Double> d1 = this.Start();
    Node<Double> d2 = seq.Start();
    int i = 0;
    for (; d1.data != null && d2.data != null; d1 = d1.next, d2 = d2.next) {
      r.Append(new Double(d1.data));
      i++;
      if (i == step) {
        r.Append(new Double(d2.data));
        i = 0;
      }
    }
    if (d1.data != null) {
      for (; d1.data != null; d1 = d1.next) {
        r.Append(new Double(d1.data));
      }
    } else if (d2.data != null) {
      for (; d2.data != null; d2 = d2.next) {
        r.Append(new Double(d2.data));
      }
    }
    return r;
  }

  //inserts at a point
  NumberSequance Inserted(NumberSequance seq, int at) {
    NumberSequance r = new NumberSequance();
    int i = 0;
    for (Node<Double> d1 = this.Start(); d1.data != null; d1 = d1.next) {
      r.Append(new Double(d1.data));
      if (i == at) {
        for (Node<Double> d2 = seq.Start(); d2.data != null; d2 = d2.next) {
          r.Append(new Double(d2.data));
        }
      }
      i++;
    }
    return r;
  }

  //gets index range
  NumberSequance Get(int from, int to) {
    NumberSequance r = new NumberSequance();
    int i = 0;
    for (Node<Double> d1 = this.Start(); d1.data != null; d1 = d1.next) {
      if (i >= from && i < to) {
        r.Append(new Double(d1.data));
      }
      i++;
    }
    return r;
  }

  //interweaves 
  NumberSequance Intergrated(NumberSequance seq) {
    NumberSequance r = new NumberSequance();
    Double s = seq.First();
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {
      r.Append(new Double(d.data));
      if (s != null) {
        r.Append(new Double(s));
        s = seq.Next();
      }
    }
    while (s != null) {
      r.Append(new Double(s));
      s = seq.Next();
    }
    return r;
  }

  //interweaves
  NumberSequance Intergrated(NumberSequance seq, int at) {
    NumberSequance r = new NumberSequance();
    Double s = seq.First();
    int i = 0;
    for (Node<Double> d = this.Start(); d.data != null; d = d.next) {

      if (i == at) {
        if (s != null) {
          r.Append(new Double(s));
          s = seq.Next();
        }
      } else {
        i++;
      }
      r.Append(new Double(d.data));
    }
    while (s != null) {
      r.Append(new Double(s));
      s = seq.Next();
    }
    return r;
  }

  NumberSequance Merge(NumberSequance seq) {
    NumberSequance r = new NumberSequance();
    r.PasteAll(this);
    r.PasteAll(seq);
    this.Clear();
    seq.Clear();
    return r;
  }

  NumberSequance Merge(NumberSequance seq[]) {
    NumberSequance r = new NumberSequance();
    r.PasteAll(this);
    for (int i = 0; i < seq.length; i++) {
      r.PasteAll(seq[i]);
      seq[i].Clear();
    }
    this.Clear();
    return r;
  }

  /*java.util.stream.Stream<Double> ToStream(NumberSequance ns) {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(ns::Get).limit(ns.length);
   return str;
   }
   
   java.util.stream.Stream<Double> ToStream() {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   return str;
   }
   
   java.util.stream.Stream<Double> ToRemoveStream() {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Remove).limit(this.length);
   return str;
   }
   
   java.util.stream.Stream<Double> ToRemoveStream(NumberSequance ns) {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(ns::Remove).limit(ns.length);
   return str;
   }
   
   void Sort(java.util.Comparator<Double> c) {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   java.util.List<Double> usl = str.sorted(c).collect(Collectors.toList());
   str.close();
   this.Clear();
   this.PasteAll(new RingList<Double>(usl));
   }
   
   ArrayList<Double> SortToList(java.util.Comparator<Double> c) {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   java.util.List<Double> usl = str.sorted(c).collect(Collectors.toList());
   str.close();
   return new ArrayList<>(usl);
   }
   
   RingList<Double> SortToRingList(java.util.Comparator<Double> c) {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   java.util.List<Double> usl = str.sorted(c).collect(Collectors.toList());
   str.close();
   return new RingList<Double>(usl);
   }
   
   NumberSequance Sorted(java.util.Comparator<Double> c) {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   java.util.List<Double> usl = str.sorted(c).collect(Collectors.toList());
   str.close();
   return new NumberSequance(usl);
   }
   
   Double[] ToDoubles() {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   java.util.List<Double> usl = str.collect(Collectors.toList());
   str.close();
   Double d[] = new Double[usl.size()];
   usl.toArray(d);
   return d;
   }
   
   Integer[] ToIntegers() {
   java.util.stream.Stream<Double> str = java.util.stream.Stream.generate(this::Get).limit(this.length);
   java.util.List<Double> usl = str.collect(Collectors.toList());
   str.close();
   Double d[] = new Double[usl.size()];
   Integer i[] = new Integer[usl.size()];
   usl.toArray(d);
   for (int n = 0; n < i.length; n++) {
   i[n] = d[n].intValue();
   }
   return i;
   }*/

  @Override
    public double f(double x) {
    if (function != null) {
      return function.f(x);
    } else {
      return x;
    }
  }

  double SumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += f(i.data);
    }
    return x;
  }

  double InvSumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += 1.0 / f(i.data);
    }
    return x;
  }

  double ProductSumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x *= f(i.data);
    }
    return x;
  }

  double InvProductSumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x /= f(i.data);
    }
    return x;
  }

  double SumLoopFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += i.data;
      x = f(x);
    }
    return x;
  }

  double InvLoopSumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {

      x = f(x);
      x += 1.0 / i.data;
    }
    return x;
  }

  double ProductLoopSumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x *= i.data;
      x = f(x);
    }
    return x;
  }

  double InvProductLoopSumFunction() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {

      x = f(x);
      x /= i.data;
    }
    return x;
  }

  double Sum() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += i.data;
    }
    return x;
  }

  double InvSum() {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += 1.0 / i.data;
    }
    return x;
  }

  double ProductSum() {
    double x = 1;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x *= i.data;
    }
    return x;
  }

  double InvProductSum() {
    double x = 1;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x /= i.data;
    }
    return x;
  }

  double Sum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += fn.f(f(i.data));
    }
    return x;
  }

  double InvSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += 1.0 / fn.f(f(i.data));
    }
    return x;
  }

  double ProductSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x *= fn.f(f(i.data));
    }
    return x;
  }

  double InvProductSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x /= fn.f(f(i.data));
    }
    return x;
  }

  double LoopSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += f(i.data);
      x = fn.f(x);
    }
    return x;
  }

  double InvLoopSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x += f(i.data);
      x = 1.0 / fn.f(x);
    }
    return x;
  }

  double LoopProductSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      x *= f(i.data);
      x = fn.f(x);
    }
    return x;
  }

  double InvLoopProductSum(OneStateFunction fn) {
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {

      x = 1.0 / fn.f(x);
      x /= f(i.data);
    }
    return x;
  }

  NumberSequance Calculate() {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data = f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateProduct() {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data *= f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateInvProduct() {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data /= f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateAll(OneStateFunction fn) {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data = fn.f(f(i.data));
      }
    }
    return this;
  }

  NumberSequance CalculateInvAll(OneStateFunction fn) {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data = f(fn.f(i.data));
      }
    }
    return this;
  }

  NumberSequance Calculate(OneStateFunction fn) {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data = fn.f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateProduct(OneStateFunction fn) {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data *= fn.f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateInvProduct(OneStateFunction fn) {
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        i.data /= fn.f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateLoop() {
    int l = Length();

    if (l > 1) {
      Double s = null;

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (s != null) {
          i.data = f(s);
        }
        s = i.data;
      }
    }
    return this;
  }

  NumberSequance CalculateLoop(OneStateFunction fn) {
    int l = Length();

    if (l > 1) {
      Double s = null;

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (s != null) {
          i.data = fn.f(s);
        }
        s = i.data;
      }
    }
    return this;
  }

  NumberSequance CalcLoop(OneStateFunction fn) {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        i.data = fn.f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalcLoop(OneStateFunction fn[]) {
    int l = Length();

    if (l > 1) {
      int j = 0;
      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        i.data = fn[j].f(i.data);
        j++;
        if (j == fn.length) {
          j = 0;
        }
      }
    }
    return this;
  }

  NumberSequance CalcLoop(OneStateFunction fn1, OneStateFunction fn2) {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
        i.data = fn1.f(i.data);
        i.next.data = fn2.f(i.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2) {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn3.f(i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn3.f2(x, y) * i.next.next.data;
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3, OneStateFunction fn4) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn4.f(fn3.f2(x, y) * i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3, TwoStateFunction fn4) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn4.f2(fn3.f2(x, y), i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, TwoStateFunction fn4) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn4.f2(x, y) + fn3.f(i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, TwoStateFunction fn4, OneStateFunction fn5) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn5.f(fn4.f2(x, y) + fn3.f(i.next.next.data));
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, TwoStateFunction fn4, TwoStateFunction fn5) {
    int l = Length();

    if (l > 2) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        double x = fn1.f(i.data) - fn2.f(i.next.data);
        double y = fn2.f(i.data) + fn1.f(i.next.data);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn5.f2(fn4.f2(x, y), fn3.f(i.next.next.data));
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, double t) {
    int l = Length();

    if (l > 1) {
      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, double t) {
    int l = Length();

    if (l > 2) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn3.f(i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3, double t) {
    int l = Length();

    if (l > 2) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn3.f2(x, y) * i.next.next.data;
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3, OneStateFunction fn4, double t) {
    int l = Length();

    if (l > 2) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn4.f(fn3.f2(x, y) * i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3, TwoStateFunction fn4, double t) {
    int l = Length();

    if (l > 2) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn4.f2(fn3.f2(x, y), i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, TwoStateFunction fn4, double t) {
    int l = Length();

    if (l > 2) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn4.f2(x, y) + fn3.f(i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, TwoStateFunction fn4, OneStateFunction fn5, double t) {
    int l = Length();

    if (l > 1) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn5.f(fn4.f2(x, y) + fn3.f(i.next.next.data));
      }
    }
    return this;
  }

  NumberSequance CalcRotaryLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3, TwoStateFunction fn4, TwoStateFunction fn5, double t) {
    int l = Length();

    if (l > 1) {

      double d = a;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        d += t;
        double x = i.data * fn1.f(d) - i.next.data * fn2.f(d);
        double y = i.data * fn2.f(d) + i.next.data * fn1.f(d);
        i.data = x;
        i.next.data = y;
        i.next.next.data = fn5.f2(fn4.f2(x, y), fn3.f(i.next.next.data));
      }
    }
    return this;
  }

  NumberSequance CalcLoop(OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3) {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
        i.data = fn1.f(i.data);
        i.next.data = fn2.f(i.next.data);
        i.next.next.data = fn3.f(i.next.next.data);
      }
    }
    return this;
  }

  NumberSequance Calc2Loop(TwoStateFunction fn1) {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
        double d = fn1.f2(i.data, i.next.data);
        i.next.data = fn1.f2(i.next.data, i.data);
        i.data = d;
      }
    }
    return this;
  }

  NumberSequance Calc2Loop(TwoStateFunction fn1, TwoStateFunction fn2) {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
        double d = fn1.f2(i.data, i.next.data);
        i.next.data = fn2.f2(i.next.data, i.data);
        i.data = d;
      }
    }
    return this;
  }

  NumberSequance Calc2Loop(TwoStateFunction fn1[]) {
    int l = Length();

    if (l > 1) {
      int j = 0;
      for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
        double d = fn1[j].f2(i.data, i.next.data);
        i.next.data = fn1[j].f2(i.next.data, i.data);
        i.data = d;
        j++;
        if (j == fn1.length) {
          j = 0;
        }
      }
    }
    return this;
  }

  NumberSequance CalcLoop() {
    int l = Length();

    if (l > 1) {

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        i.data = f(i.data);
      }
    }
    return this;
  }

  NumberSequance CalculateLoopSum(OneStateFunction fn) {
    int l = Length();

    if (l > 1) {
      Double s = null;

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (s != null) {
          i.data += fn.f(s);
        }
        s = i.data;
      }
    }
    return this;
  }

  NumberSequance CalculateLoopAll(OneStateFunction fn) {
    int l = Length();

    if (l > 1) {
      Double s = null;

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (s != null) {
          i.data = fn.f(f(s));
        }
        s = i.data;
      }
    }
    return this;
  }

  NumberSequance CalculateLoopProduct(OneStateFunction fn) {
    int l = Length();

    if (l > 1) {
      Double s = null;

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (s != null) {
          i.data *= fn.f(s);
        }
        s = i.data;
      }
    }
    return this;
  }

  NumberSequance CalculateLoopInvProduct(OneStateFunction fn) {
    int l = Length();

    if (l > 1) {
      Double s = null;

      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (s != null) {
          i.data /= fn.f(s);
        }
        s = i.data;
      }
    }
    return this;
  }

  NumberSequance Operation(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn.f2(i.data, d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Operation(OneStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      copy.Append(fn.f(i.data));
    }

    return copy;
  }

  NumberSequance OperationAll(OneStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      copy.Append(f(fn.f(i.data)));
    }

    return copy;
  }

  NumberSequance OperationInvAll(OneStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      copy.Append(fn.f(f(i.data)));
    }

    return copy;
  }

  NumberSequance Operation(TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
      copy.Append(fn.f2(i.data, i.next.data));
    }

    return copy;
  }

  NumberSequance Operation(ThreeStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
      copy.Append(fn.f3(i.data, i.next.data, i.next.next.data));
    }

    return copy;
  }

  NumberSequance Operation(FourStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null && i.next.next.next.data != null; i = i.next.next.next.next) {
      copy.Append(fn.f4(new Vector(i.data, i.next.data, i.next.next.data, i.next.next.next.data)));
    }
    return copy;
  }

  NumberSequance OperationAll(TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
      copy.Append(f(fn.f2(i.data, i.next.data)));
    }

    return copy;
  }

  NumberSequance OperationAll(ThreeStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
      copy.Append(f(fn.f3(i.data, i.next.data, i.next.next.data)));
    }

    return copy;
  }

  NumberSequance OperationAll(FourStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null && i.next.next.next.data != null; i = i.next.next.next.next) {
      copy.Append(f(fn.f4(new Vector(i.data, i.next.data, i.next.next.data, i.next.next.next.data))));
    }
    return copy;
  }

  NumberSequance OperationEveryOther(OneStateFunction fn1, OneStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    boolean alt = true;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      if (alt) {
        copy.Append(fn1.f(i.data));
      } else {
        copy.Append(fn2.f(i.data));
      }
      alt = !alt;
    }
    return copy;
  }

  NumberSequance OperationEveryOther(OneStateFunction fn[]) {
    NumberSequance copy = new NumberSequance();
    int alt = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      copy.Append(fn[alt].f(i.data));
      alt++;
      if (alt == fn.length) {
        alt = 0;
      }
    }
    return copy;
  }

  NumberSequance Operation(NumberSequance seq, TwoStateFunction fn1, TwoStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        boolean alt = true;
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          if (alt) {
            copy.Append(fn1.f2(i.data, d));
          } else {
            copy.Append(fn2.f2(i.data, d));
          }
          alt = !alt;
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Operation(NumberSequance seq, TwoStateFunction fn[]) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        int alt = 0;
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn[alt].f2(i.data, d));
          alt++;
          if (alt == fn.length) {
            alt = 0;
          }
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Addition(NumberSequance seq) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data + d);
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Addition(NumberSequance seq, OneStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data + fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Addition(NumberSequance seq, OneStateFunction fn1, OneStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn1.f(i.data) + fn2.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Addition(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data + fn.f2(i.data, d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance AdditionInv(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(d + fn.f2(d, i.data));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance AdditionAll(NumberSequance seq, OneStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(f(i.data) + fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance AdditionInvAll(NumberSequance seq, OneStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn.f(i.data) + f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Subtract(NumberSequance seq) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data - d);
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Subtract(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data - fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Subtract(NumberSequance seq, OneStateFunction fn1, OneStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn1.f(i.data) - fn2.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Subtract(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data - fn.f2(i.data, d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance SubtractInv(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(d - fn.f2(d, i.data));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance SubtractAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(f(i.data) - fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance SubtractInvAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn.f(i.data) - f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Multiply(NumberSequance seq) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data * d);
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Multiply(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data * fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Multiply(NumberSequance seq, OneStateFunction fn1, OneStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn1.f(i.data) * fn2.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Multiply(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data * fn.f2(i.data, d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance MultipyInv(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(d * fn.f2(d, i.data));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance MultiplyAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(f(i.data) * fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance MultiplyInvAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn.f(i.data) * f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Divide(NumberSequance seq) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data / d);
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Divide(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data / fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Divide(NumberSequance seq, OneStateFunction fn1, OneStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn1.f(i.data) / fn2.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Divide(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(i.data / fn.f2(i.data, d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance DivideInv(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(d / fn.f2(d, i.data));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance DivideAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(f(i.data) / fn.f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance DivideInvAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(fn.f(i.data) / f(d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Power(NumberSequance seq) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(i.data, d));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Power(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(i.data, fn.f(d)));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Power(NumberSequance seq, OneStateFunction fn1, OneStateFunction fn2) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(fn1.f(i.data), fn2.f(d)));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance Power(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(i.data, fn.f2(i.data, d)));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance PowerInv(NumberSequance seq, TwoStateFunction fn) {
    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(d, fn.f2(d, i.data)));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance PowerAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(f(i.data), fn.f(d)));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  NumberSequance PowerInvAll(NumberSequance seq, OneStateFunction fn) {

    NumberSequance copy = new NumberSequance();
    if (seq != null) {
      Double d = seq.First();
      if (d != null) {
        for (Node<Double> i = Start(); i.data != null; i = i.next) {
          copy.Append(Math.pow(fn.f(i.data), f(d)));
          d = seq.Next();
          if (d == null) {
            d = seq.First();
          }
        }
      }
    }
    return copy;
  }

  double OperateSum() {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data += f(i.data);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateSumLoop(double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x = i.data += f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopSum(double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data += f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopProduct(double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data = i.data * f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopInvProduct(double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data = i.data / f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateInvProduct() {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += (i.data /= f(i.data));
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateSum(OneStateFunction fn) {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data += fn.f(i.data);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopSum(OneStateFunction fn) {
    Double x = null;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {
        if (x != null) {
          i.data += fn.f(x);
        }
        x = i.data;
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopProduct(OneStateFunction fn) {
    Double x = null;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        if (x != null) {
          i.data = i.data * fn.f(x);
        }
        x = i.data;
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopInvProduct(OneStateFunction fn) {
    Double x = null;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        if (x != null) {
          i.data = i.data / fn.f(x);
        }
        x = i.data;
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopSum() {
    Double x = null;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        if (x != null) {
          i.data += f(x);
        }
        x = i.data;
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopProduct() {
    Double x = null;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        if (x != null) {
          i.data = i.data * f(x);
        }
        x = i.data;
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopInvProduct() {
    Double x = null;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        if (x != null) {
          i.data = i.data / f(x);
        }
        x = i.data;
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateProduct() {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += (i.data *= f(i.data));
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopSum(OneStateFunction fn, double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data += fn.f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopProduct(OneStateFunction fn, double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data = i.data * fn.f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateLoopInvProduct(OneStateFunction fn, double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data = i.data / fn.f(x);
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateProduct(OneStateFunction fn) {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += (i.data *= fn.f(i.data));
      }
      return x;
    } else {
      return x;
    }
  }

  double OperateInvProduct(OneStateFunction fn) {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += (i.data /= fn.f(i.data));
      }
      return x;
    } else {
      return x;
    }
  }

  double Average() {
    double x = 0;
    int l = Length();

    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += i.data;
      }
      return x / l;
    } else {
      return x;
    }
  }

  double AverageFunction() {
    double x = 0;
    int l = Length();
    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += f(i.data);
      }
      return x / l;
    } else {
      return x;
    }
  }

  double AverageFunctionAll(OneStateFunction fn) {
    double x = 0;
    int l = Length();
    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += fn.f(f(i.data));
      }
      return x / l;
    } else {
      return x;
    }
  }

  double AverageInvFunctionAll(OneStateFunction fn) {
    double x = 0;
    int l = Length();
    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += f(fn.f(i.data));
      }
      return x / l;
    } else {
      return x;
    }
  }

  double AverageFunction(OneStateFunction fn) {
    double x = 0;
    int l = Length();
    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += fn.f(i.data);
      }
      return x / l;
    } else {
      return x;
    }
  }

  double AverageFunctionResult(OneStateFunction fn) {
    double x = 0;
    int l = Length();
    if (l > 0) {
      for (Node<Double> i = Start(); i.data != null; i = i.next) {

        x += f(i.data);
      }
      return fn.f(x) / l;
    } else {
      return fn.f(x);
    }
  }

  RingList<Vector> ToVectors2D() {
    RingList<Vector> vectors = new RingList();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
      Vector v = new Vector(i.data, i.next.data);
      vectors.Append(v);
    }
    return vectors;
  }

  RingList<Vector> ToVectors3D() {
    RingList<Vector> vectors = new RingList();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
      Vector v = new Vector(i.data, i.next.data, i.next.next.data);
      vectors.Append(v);
    }
    return vectors;
  }

  RingList<Point> ToPoints2D() {
    RingList<Point> vectors = new RingList();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
      Point v = new Point(i.data, i.next.data);
      vectors.Append(v);
    }
    return vectors;
  }

  RingList<Point> ToPoints3D() {
    RingList<Point> vectors = new RingList();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
      Point v = new Point(i.data, i.next.data, i.next.next.data);
      vectors.Append(v);
    }
    return vectors;
  }

  RingList<Point> ToRadialPoints3D() {
    RingList<Point> vectors = new RingList();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null && i.next.next.next.data != null; i = i.next.next.next.next) {
      Point v = new Point(i.data, i.next.data, i.next.next.data, i.next.next.next.data);
      vectors.Append(v);
    }
    return vectors;
  }

  RingList<Point> ToRadialPoints2D() {
    RingList<Point> vectors = new RingList();
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
      Point v = new Point(i.data, i.next.data);
      v.r = i.next.next.data;
      vectors.Append(v);
    }
    return vectors;
  }

  NumberSequance GenerateLoop(int number, double seed) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x = f(x);
      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateSumLoop(NumberSequance seq, double seed) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = i.data + f(x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoop(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = fn.f2(i.data, x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopInv(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = fn.f2(x, i.data);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopInvProduct(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x /= fn.f2(i.data, x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopSum(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x += fn.f2(i.data, x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopSumInv(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x += fn.f2(x, i.data);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopProduct(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x *= fn.f2(i.data, x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopAll(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = fn.f2(f(i.data), x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateAllLoop(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = f(fn.f2(i.data, x));

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateAllLoopInv(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = f(fn.f2(x, i.data));

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopInvAll(NumberSequance seq, double seed, TwoStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = fn.f2(i.data, f(x));

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopSum(int number, double seed) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x += f(x);
      this.Append(x);
    }
    return this;
  }
  //

  NumberSequance GenerateProductLoop(NumberSequance seq, double seed) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = i.data * f(x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopProduct(int number, double seed) {
    double x = seed;
    int l = Length();

    if (l > 0) {
      for (int i = 0; i < number; i++) {

        x *= f(x);

        this.Append(x);
      }
    }
    return this;
  }

  NumberSequance GenerateProductInvLoop(NumberSequance seq, double seed) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = i.data / f(x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopInvProduct(int number, double seed) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x /= f(x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed, OneStateFunction fn) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x = fn.f(x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed, OneStateFunction fn1[]) {
    double x = seed;
    int j = 0;
    for (int i = 0; i < number; i++) {

      x = fn1[j].f(x);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, double seed3, ThreeStateFunction fn1[]) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    int j = 0;
    for (int i = 0; i < number; i++) {

      double x1 = fn1[j].f3(x, y, z);
      double y1 = fn1[j].f3(y, z, x);
      double z1 = fn1[j].f3(z, x, y);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(x1);
      this.Append(y1);
      this.Append(z1);
      x = x1;
      y = y1;
      z = z1;
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, double seed3, ThreeStateFunction fn1[], OneStateFunction f1) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    int j = 0;
    for (int i = 0; i < number; i++) {
      double x1 = fn1[j].f3(x, y, z) * f1.f(x);
      double y1 = fn1[j].f3(y, z, x) * f1.f(y);
      double z1 = fn1[j].f3(z, x, y) * f1.f(z);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(x1);
      this.Append(y1);
      this.Append(z1);
      x = x1;
      y = y1;
      z = z1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, OneStateFunction f1, OneStateFunction f2) {
    double x = seed1;
    double y = seed2;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);

      this.Append(x1);
      this.Append(y1);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);
      z = f3.f(z);

      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, TwoStateFunction f3) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);
      z = f3.f2(x, y) * z;

      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, TwoStateFunction f4) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);
      z = f3.f(z) + f4.f2(x, y);
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, TwoStateFunction f3, TwoStateFunction f4) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);
      z = f3.f2(f4.f2(x, y), z);
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, TwoStateFunction f4, OneStateFunction f5) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);
      z = f5.f(f3.f(z) + f4.f2(x, y));
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, TwoStateFunction f4, TwoStateFunction f5) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    for (int i = 0; i < number; i++) {

      double x1 = f1.f(x) - f2.f(y);
      double y1 = f1.f(y) + f2.f(x);
      z = f5.f2(f4.f2(x, y), f3.f(z));
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, OneStateFunction f1, OneStateFunction f2, double t) {
    double x = seed1;
    double y = seed2;
    double d = a;
    for (int i = 0; i < number; i++) {
      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);

      this.Append(x1);
      this.Append(y1);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, double t) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    double d = a;
    for (int i = 0; i < number; i++) {

      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);
      z = f3.f(z);

      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, TwoStateFunction f3, double t) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    double d = a;
    for (int i = 0; i < number; i++) {

      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);
      z = f3.f2(x1, y1) * z;

      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, TwoStateFunction f4, double t) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    double d = a;
    for (int i = 0; i < number; i++) {

      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);
      z = f3.f(z) + f4.f2(x1, y1);
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, TwoStateFunction f3, TwoStateFunction f4, double t) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    double d = a;
    for (int i = 0; i < number; i++) {

      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);
      z = f3.f2(f4.f2(x1, y1), z);
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, TwoStateFunction f4, OneStateFunction f5, double t) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    double d = a;
    for (int i = 0; i < number; i++) {

      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);
      z = f5.f(f3.f(z) + f4.f2(x1, y1));
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateRotaryLoop(int number, double seed1, double seed2, double seed3, OneStateFunction f1, OneStateFunction f2, OneStateFunction f3, TwoStateFunction f4, TwoStateFunction f5, double t) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    double d = a;
    for (int i = 0; i < number; i++) {

      d += t;
      double x1 = x * f1.f(d) - y * f2.f(d);
      double y1 = y * f1.f(d) + x * f2.f(d);
      z = f5.f2(f4.f2(x1, y1), f3.f(z));
      this.Append(x1);
      this.Append(y1);
      this.Append(z);
      x = x1;
      y = y1;
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, OneStateFunction fn1[]) {
    double x = seed1;
    double z = seed2;
    int j = 0;
    for (int i = 0; i < number; i += 2) {

      x = fn1[j].f(x);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(x);
      z = fn1[j].f(z);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(z);
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, double seed3, OneStateFunction fn1[]) {
    double x = seed1;
    double y = seed2;
    double z = seed3;
    int j = 0;
    for (int i = 0; i < number; i += 3) {

      x = fn1[j].f(x);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(x);
      y = fn1[j].f(y);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(y);
      z = fn1[j].f(z);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
      this.Append(z);
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, TwoStateFunction fn1[]) {
    double x = seed1;
    int j = 0;
    double z = seed2;
    for (int i = 0; i < number; i += 2) {

      x = fn1[j].f2(x, z);
      this.Append(x);
      z = fn1[j].f2(z, x);
      this.Append(z);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, TwoStateFunction fn1[], TwoStateFunction fn2[]) {
    double x = seed1;
    int j1 = 0;
    int j2 = 0;
    double z = seed2;
    for (int i = 0; i < number; i += 2) {

      x = fn1[j1].f2(x, z);
      this.Append(x);
      z = fn2[j2].f2(z, x);
      this.Append(z);
      j1++;
      if (j1 == fn1.length) {
        j1 = 0;
      }
      j2++;
      if (j2 == fn2.length) {
        j2 = 0;
      }
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, double seed3, TwoStateFunction fn1[]) {
    double x = seed1;
    double y = seed2;
    int j = 0;
    double z = seed3;
    for (int i = 0; i < number; i += 3) {

      x = fn1[j].f2(z, x);
      this.Append(x);
      y = fn1[j].f2(x, y);
      this.Append(y);
      z = fn1[j].f2(y, z);
      this.Append(z);
      j++;
      if (j == fn1.length) {
        j = 0;
      }
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed1, double seed2, double seed3, TwoStateFunction fn1[], TwoStateFunction fn2[], TwoStateFunction fn3[]) {
    double x = seed1;
    double y = seed2;
    int j1 = 0;
    int j2 = 0;
    int j3 = 0;
    double z = seed3;
    for (int i = 0; i < number; i += 2) {

      x = fn1[j1].f2(z, x);
      this.Append(x);
      y = fn2[j2].f2(x, y);
      this.Append(y);
      z = fn3[j3].f2(y, z);
      this.Append(z);
      j1++;
      if (j1 == fn1.length) {
        j1 = 0;
      }
      j2++;
      if (j2 == fn2.length) {
        j2 = 0;
      }
      j3++;
      if (j3 == fn3.length) {
        j3 = 0;
      }
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed, OneStateFunction fn1, OneStateFunction fn2) {
    double x = seed;
    for (int i = 0; i < number; i += 2) {

      x = fn1.f(x);
      this.Append(x);
      x = fn2.f(x);
      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed, OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x = fn3.f2(fn1.f(x), fn2.f(x));
      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopCombinations(int number, double seed, OneStateFunction fn1, OneStateFunction fn2, TwoStateFunction fn3) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      double y = fn3.f2(fn1.f(x), fn1.f(x));
      this.Append(x);
      y = fn3.f2(fn1.f(x), fn2.f(x));
      this.Append(x);
      y = fn3.f2(fn2.f(x), fn2.f(x));
      this.Append(x);
      y = fn3.f2(fn2.f(x), fn1.f(x));
      this.Append(x);
      x = y;
    }
    return this;
  }

  NumberSequance GenerateLoop(int number, double seed, OneStateFunction fn1, OneStateFunction fn2, OneStateFunction fn3) {
    double x = seed;
    for (int i = 0; i < number; i += 3) {

      x = fn1.f(x);
      this.Append(x);
      x = fn2.f(x);
      this.Append(x);
      x = fn3.f(x);
      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateSumLoop(NumberSequance seq, double seed, OneStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = i.data + fn.f(x);
      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopSum(int number, double seed, OneStateFunction fn) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x += fn.f(x);
      this.Append(x);
    }
    return this;
  }
  //

  NumberSequance GenerateProductLoop(NumberSequance seq, double seed, OneStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = i.data * fn.f(x);

      this.Append(x);
    }
    return this;
  }

  NumberSequance GenerateLoopProduct(int number, double seed, OneStateFunction fn) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x *= fn.f(x);

      this.Append(x);
    }

    return this;
  }

  NumberSequance GenerateProductInvLoop(NumberSequance seq, double seed, OneStateFunction fn) {
    double x = seed;
    for (Node<Double> i = seq.Start(); i.data != null; i = i.next) {

      x = i.data / fn.f(x);
      this.Append(x);
    }

    return this;
  }

  NumberSequance GenerateLoopInvProduct(int number, double seed, OneStateFunction fn) {
    double x = seed;
    for (int i = 0; i < number; i++) {

      x /= fn.f(x);
      this.Append(x);
    }
    return this;
  }

  NumberSequance Generate(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    while (n >= from && n < to) {
      Double res = fn.f(n);
      n += inc;
      this.Append(res);
    }
    return this;
  }

  NumberSequance GenerateSum(int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = 0;
    while (n >= from && n < to) {
      res += f(n);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateProduct(int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = 0;
    while (n >= from && n < to) {
      res *= f(n);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateLoop(int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = n;
    while (n >= from && n < to) {
      res = f(res);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateProductLoop(int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = n;
    while (n >= from && n < to) {
      res *= f(res);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateInvProductLoop(int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = n;
    while (n >= from && n < to) {
      res /= f(res);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateInvProduct(int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = 0;
    while (n >= from && n < to) {
      res /= f(n);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateSum(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = 0;
    while (n >= from && n < to) {
      res += fn.f(n);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateProduct(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = 0;
    while (n >= from && n < to) {
      res *= fn.f(n);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateLoop(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = n;
    while (n >= from && n < to) {
      res = fn.f(res);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateProductLoop(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = n;
    while (n >= from && n < to) {
      res *= fn.f(res);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateInvProductLoop(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = n;
    while (n >= from && n < to) {
      res /= fn.f(res);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  NumberSequance GenerateInvProduct(OneStateFunction fn, int from, int to, int inc) {
    a = from;
    b = to;
    n = from;
    step = inc;
    if ((from > to && inc > 0) || (from < to && inc < 0)) {
      int s = from;
      from = to;
      to = s;
    }
    double res = 0;
    while (n >= from && n < to) {
      res /= fn.f(n);
      n += inc;
      this.Append(new Double(res));
    }
    return this;
  }

  String ToString(int lineFeed) {
    StringBuilder str = new StringBuilder();
    int j = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      str.append(j);
      str.append("=" + (i.data));
      j++;

      if (j % lineFeed == 0) {
        str.append("\n");
      } else {
        str.append("\t");
      }
    }
    return str.toString();
  }
  String ToStringSum(int lineFeed) {
    StringBuilder str = new StringBuilder();
    int j = 0;
    double x = 0;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      str.append(j);
      str.append("=" + (i.data.floatValue()));
      x += i.data;
      j++;

      if (j % lineFeed == 0) {
        str.append("\n");
      } else {
        str.append("\t");
      }
    }
    str.append("\nSum= ");
    str.append(x);
    str.append("\n");
    return str.toString();
  }

  String ToStringXY() {
    StringBuilder str = new StringBuilder();
    int j = 0;
    double x = 0;
    double max = 0;
    for (Node<Double> i = Start(); i.data != null && i.next.data != null; i = i.next.next) {
      str.append(j);
      str.append(" x").append(i.data.floatValue());
      str.append(" y").append(i.next.data.floatValue());
      x = i.data;
      if (Math.abs(x) > max) {
        max = Math.abs(x);
      }
      x = i.next.data;
      if (Math.abs(x) > max) {
        max = Math.abs(x);
      }
      j++;

      if (j % 2 == 0) {
        str.append("\n");
      } else {
        str.append("\t");
      }
    }
    str.append("\nRadius= ");
    str.append(x);
    str.append("\n");
    return str.toString();
  }

  String ToStringXYZ() {
    StringBuilder str = new StringBuilder();
    int j = 0;
    double x = 0;
    double max = 0;
    for (Node<Double> i = Start(); i.data != null && i.next.data != null && i.next.next.data != null; i = i.next.next.next) {
      str.append(j);
      str.append(" x").append(i.data.floatValue());
      str.append(" y").append(i.next.data.floatValue());
      str.append(" z").append(i.next.next.data.floatValue());
      x = i.data;
      if (Math.abs(x) > max) {
        max = Math.abs(x);
      }
      x = i.next.data;
      if (Math.abs(x) > max) {
        max = Math.abs(x);
      }
      x = i.next.next.data;
      if (Math.abs(x) > max) {
        max = Math.abs(x);
      }
      j++;

      if (j % 3 == 0) {
        str.append("\n");
      } else {
        str.append("\t");
      }
    }
    str.append("\nRadius= ");
    str.append(x);
    str.append("\n");
    return str.toString();
  }

  Limit CalcRange() {
    Limit lim = new Limit();
    lim.from = Double.MAX_VALUE;
    lim.to = Double.MIN_VALUE;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      if (i.data < lim.to) {
        lim.to = i.data;
      } else if (i.data > lim.from) {
        lim.from = i.data;
      }
    }
    return lim;
  }

  Limit CalcRange(OneStateFunction fn) {
    Limit lim = new Limit();
    lim.from = Double.MAX_VALUE;
    lim.to = Double.MIN_VALUE;
    for (Node<Double> i = Start(); i.data != null; i = i.next) {
      double data = fn.f(i.data);
      if (data < lim.to) {
        lim.to = data;
      } else if (data > lim.from) {
        lim.from = data;
      }
    }
    return lim;
  }

  Limit CalcLimit() {
    Limit lim = new Limit();
    lim.from = Double.MIN_VALUE;
    double past = 0;
    Node<Double> i = Start();
    if (i.data != null) {
      double data = past = i.data;

      lim.from = data;
      i = i.next;
      for (; i.data != null; i = i.next) {
        data = i.data;
        if (data < past) {
          past = data;
        } else if (data > past) {
          past = data;
        }
      }
      lim.to = past;
    }
    return lim;
  }

  Limit CalcLimit(OneStateFunction fn) {
    Limit lim = new Limit();
    lim.from = Double.MIN_VALUE;
    double past = 0;
    Node<Double> i = Start();
    if (i.data != null) {
      double data = past = fn.f(i.data);

      lim.from = data;
      i = i.next;
      for (; i.data != null; i = i.next) {
        data = fn.f(i.data);
        if (data < past) {
          past = data;
        } else if (data > past) {
          past = data;
        }
      }
      lim.to = past;
    }
    return lim;
  }

  Limit CalcAverageLimit() {
    Limit lim = new Limit();
    lim.from = Double.MIN_VALUE;
    lim.to = Double.MAX_VALUE;
    Node<Double> i = Start();
    double past = 0;
    if (i.data != null) {
      double data = past = i.data;

      lim.from = data;
      i = i.next;
      int j = 1;
      for (; i.data != null; i = i.next) {
        if (i.data * j < past) {
          past += i.data;
          j++;
        } else if (i.data * j > past) {
          past += i.data;
          j++;
        }
      }
      lim.to = past / this.Length();
    }
    return lim;
  }

  Limit CalcAverageLimit(OneStateFunction fn) {
    Limit lim = new Limit();
    lim.from = Double.MIN_VALUE;
    lim.to = Double.MAX_VALUE;
    Node<Double> i = Start();
    double past = 0;
    if (i.data != null) {
      double data = past = i.data;

      lim.from = data;
      i = i.next;
      int j = 1;
      for (; i.data != null; i = i.next) {
        data = i.data;
        if (data * j < past) {
          past += data;
          j++;
        } else if (data * j > past) {
          past += data;
          j++;
        }
      }
      lim.to = past / this.Length();
    }
    return lim;
  }
}

class Limit {

  Double from;
  Double to;
}

class IntegralOperation extends RingList<NumberSequance> {

  NumberSequance variables = new NumberSequance();//subject
  NumberSequance results = new NumberSequance();//subject
  NumberSequance changes = new NumberSequance();//subject
  double start;//counter record
  double end;//counter record
  double a;//counters
  double b;//counters
  double c;//constant
  RingList<OneStateFunction> operations = new RingList();
  RingList<OneStateFunction> inverseOperations = new RingList();
  OneStateFunction iterate;// = (x) -> x;
  OneStateFunction resultant;// = (x) -> x;
  OneStateFunction intergrate=new OneStateFunction() {
    double f(double x) {
      return 0.5*x*x;
    }
  };// = (x) -> 0.5 * x * x;
  OneStateFunction invIntergrate;// = (x) -> 1.0 / x;
  //
  IntegralOperation recursive;//functions only
  NumberSequance record;
  boolean recording = false;

  IntegralOperation(double a, double b, double c) {
    this.a = a;
    this.b = b;
    this.c = c;
    start = a;
    end = b;
  }

  IntegralOperation(double a, double b, double c, OneStateFunction fn) {
    this.a = a;
    this.b = b;
    this.c = c;
    operations.Append(fn);
    start = a;
    end = b;
  }

  IntegralOperation(double a, double b, double c, OneStateFunction fn, OneStateFunction inv) {
    this.a = a;
    this.b = b;
    this.c = c;
    operations.Append(fn);
    inverseOperations.Append(inv);
    start = a;
    end = b;
  }

  IntegralOperation(double a, double b, double c, OneStateFunction fn, IntegralOperation recure) {

    this.a = a;
    this.b = b;

    this.c = c;
    operations.Append(fn);
    recursive = recure;
    start = a;
    end = b;
  }

  IntegralOperation(double a, double b, double c, OneStateFunction fn, OneStateFunction inv, IntegralOperation recure) {
    this.a = a;
    this.b = b;
    this.c = c;
    operations.Append(fn);
    inverseOperations.Append(inv);
    recursive = recure;
    start = a;
    end = b;
  }

  void ResetState(double from, double to) {
    start = from;
    end = to;
    a = start;
    b = end;
  }

  void ResetState() {
    a = start;
    b = end;
  }

  void Reset(Double[] x) {
    variables = new NumberSequance();
    variables.Append(x);
    a = start;
    b = end;
    results = new NumberSequance();
    changes = new NumberSequance();
  }

  void Reset(RingList<Double> x) {
    variables = new NumberSequance();
    variables.Append(x);
    a = start;
    b = end;
    results = new NumberSequance();
    changes = new NumberSequance();
  }

  void Record(NumberSequance seq) {
    this.record = seq;
  }

  IntegralOperation IterateUsing(OneStateFunction fn) {
    iterate=fn;
    return this;
  }
  IntegralOperation ResultantUsing(OneStateFunction fn) {
    iterate=fn;
    return this;
  }
  IntegralOperation ItergrateUsing(OneStateFunction fn) {
    intergrate=fn;
    return this;
  }
  IntegralOperation InvItergrateUsing(OneStateFunction fn) {
    invIntergrate=fn;
    return this;
  }
  IntegralOperation OperateUsing(OneStateFunction fn) {
    operations.Append(fn);
    return this;
  }

  IntegralOperation ReOperateUsing(OneStateFunction fn) {
    operations = new RingList();
    operations.Append(fn);
    return this;
  }

  IntegralOperation ReInvertUsing(OneStateFunction fn) {
    inverseOperations = new RingList();
    inverseOperations.Append(fn);
    return this;
  }

  IntegralOperation InvertUsing(OneStateFunction fn) {
    inverseOperations.Append(fn);
    return this;
  }

  IntegralOperation Recure(IntegralOperation s) {
    if (recursive != null) {
      return recursive.Recure(s);
    } else {
      recursive = s;
    }
    return recursive;
  }

  void MergeVariables(NumberSequance seq) {
    variables = variables.Merge(seq);
  }

  void SetVariables(NumberSequance seq) {
    variables = seq;
  }

  void IntergrateVariables(NumberSequance seq) {
    variables = variables.Intergrated(seq);
  }

  void IntergrateVariables(NumberSequance seq, int at) {
    variables = variables.Intergrated(seq, at);
  }

  void CombineVariables(NumberSequance seq) {
    variables = variables.Combined(seq);
  }

  void CombineVariables(NumberSequance seq1, NumberSequance seq2) {
    variables = variables.Combined(seq1, 1);
    variables = variables.Combined(seq2, 2);
  }

  void CombineVariables(NumberSequance seq1, NumberSequance seq2, NumberSequance seq3) {
    variables = variables.Combined(seq1, 1);
    variables = variables.Combined(seq2, 2);
    variables = variables.Combined(seq3, 3);
  }

  void CombineVariables(NumberSequance seq[]) {
    for (int i = 0; i < seq.length; i++) {
      variables = variables.Combined(seq[i], i + 1);
    }
  }

  void AddVariable(Double v) {
    variables.Append(v);
  }

  void AddVariables(Double v[]) {
    variables.Append(v);
  }

  void AddVariables(RingList<Double> v) {
    variables.Append(v);
  }

  NumberSequance Variables() {
    return variables;
  }

  double Area() {
    OneStateFunction operate = operations.First();
    double oa = 0;
    while (operate != null) {
      if (recursive == null) {
        double area = 0.5 * (b - a) * (operate.f(a) + operate.f(b));
        oa += area;
      } else {
        double area = 0.5 * (b - a) * (Recure(operate, a) + Recure(operate, b));
        oa += area;
      }
      operate = operations.Next();
    }
    return oa;
  }

  double Area(double ax, double bx) {
    OneStateFunction operate = operations.First();
    double oa = 0;
    while (operate != null) {
      if (recursive == null) {
        double area = 0.5 * (bx - ax) * (operate.f(ax) + operate.f(bx));
        oa += area;
      } else {
        double area = 0.5 * (bx - ax) * (Recure(operate, ax) + Recure(operate, bx));
        oa += area;
      }
      operate = operations.Next();
    }
    return oa;
  }

  double Recure(OneStateFunction fn, double r) {
    if (recursive != null && recursive.operations.Length() > 0) {
      OneStateFunction operate = recursive.operations.First();

      double o = fn.f(r);
      while (operate != null) {
        o = Recure(operate, o);
        operate = recursive.operations.Next();
      }
      return o;
    } else {
      return r;
    }
  }

  double Recure(double x) {
    if (operations.Length() > 0) {
      OneStateFunction operate = operations.First();

      double o = operate.f(x);
      while (operate != null) {
        o = Recure(operate, o);
        operate = operations.Next();
      }
      return o;
    } else {
      return x;
    }
  }

  double RecursiveArea(double ax, double bx) {
    return 0.5 * (bx - ax) * (Recure(ax) + Recure(bx));
  }

  double RecursiveArea() {
    return 0.5 * (b - a) * (Recure(a) + Recure(b));
  }

  double InvRecure(OneStateFunction fn, double r) {
    if (recursive != null && recursive.inverseOperations.Length() > 0) {
      OneStateFunction operate = recursive.inverseOperations.Last();

      double o = r;
      while (operate != null) {
        o = InvRecure(operate, o);
        operate = recursive.inverseOperations.Prev();
      }
      return fn.f(o);
    } else {
      return r;
    }
  }

  double InvRecure(double x) {
    if (inverseOperations.Length() > 0) {
      OneStateFunction operate = inverseOperations.Last();

      double o = x;
      while (operate != null) {
        o = InvRecure(operate, o);
        operate = inverseOperations.Prev();
      }
      return operate.f(o);
    } else {
      return x;
    }
  }

  //Random Chaos Function
  double Scramble(OneStateFunction fn, double r) {
    if (recursive != null && recursive.operations.Length() > 0) {
      OneStateFunction operate = recursive.operations.First();
      OneStateFunction inverse = recursive.inverseOperations.Last();

      double o = fn.f(r);
      while (operate != null) {
        int i = arcRandomInt(4);
        switch (i) {
        case 0: 
          {
            o = Scramble(operate, o);
            break;
          }
        case 1: 
          {
            o = Scramble(inverse, o);
            break;
          }
        case 2: 
          {
            o = InvScramble(operate, o);
            break;
          }
        case 3: 
          {
            o = InvScramble(inverse, o);
            break;
          }
        }
        operate = operations.Next();
        inverse = inverseOperations.Prev();
      }
      return o;
    } else {
      return r;
    }
  }

  //Random Chaos Function
  double InvScramble(OneStateFunction fn, double r) {
    if (recursive != null && recursive.operations.Length() > 0) {
      OneStateFunction operate = recursive.operations.First();
      OneStateFunction inverse = recursive.inverseOperations.Last();

      double o = r;
      while (operate != null) {
        int i = arcRandomInt(4);
        switch (i) {
        case 0: 
          {
            o = Scramble(operate, o);
            break;
          }
        case 1: 
          {
            o = Scramble(inverse, o);
            break;
          }
        case 2: 
          {
            o = InvScramble(operate, o);
            break;
          }
        case 3: 
          {
            o = InvScramble(inverse, o);
            break;
          }
        }
        operate = operations.Next();
        inverse = inverseOperations.Prev();
      }
      return fn.f(o);
    } else {
      return r;
    }
  }

  //Random Chaos Function
  double Scramble(double x) {
    if (operations.Length() > 0) {
      OneStateFunction operate = operations.First();
      OneStateFunction inverse = inverseOperations.Last();

      double o = operate.f(x);
      while (operate != null) {
        int i = arcRandomInt(4);
        switch (i) {
        case 0: 
          {
            o = Scramble(operate, o);
            break;
          }
        case 1: 
          {
            o = Scramble(inverse, o);
            break;
          }
        case 2: 
          {
            o = InvScramble(operate, o);
            break;
          }
        case 3: 
          {
            o = InvScramble(inverse, o);
            break;
          }
        }
        operate = operations.Next();
        inverse = inverseOperations.Prev();
      }
      return o;
    } else {
      return x;
    }
  }

  boolean HasVariables() {
    return variables.Length() > 0;
  }

  boolean HasInverse() {
    return inverseOperations.Length() > 0;
  }

  boolean CanOperate() {
    return operations.Length() > 0;
  }

  int Inversions() {
    return inverseOperations.Length();
  }

  int Operators() {
    return operations.Length();
  }

  boolean CanRecord() {
    return record != null;
  }

  NumberSequance Recorded() {
    return record;
  }

  NumberSequance Results() {
    return results;
  }

  NumberSequance Changes() {
    return changes;
  }

  void StoreRecorded() {
    Append(record);
    record = new NumberSequance();
  }

  void StoreResults() {
    Append(results.Clone());
  }

  void StoreChanges() {
    Append(changes.Clone());
  }

  void StoreVariables() {
    Append(variables.Clone());
  }

  void ReStoreToVariables() {
    this.variables = this.Remove();
  }

  void ReStoreToResults() {
    this.results = this.Remove();
  }

  void ReStoreToChanges() {
    this.changes = Remove();
  }

  void ReStoreToRecord() {
    this.record = Remove();
  }

  NumberSequance CopyResults() {
    return (results.Clone());
  }

  NumberSequance CopyChanges() {
    return (changes.Clone());
  }

  NumberSequance CopyVariables() {
    return (variables.Clone());
  }

  NumberSequance CopyRecord() {
    return (record.Clone());
  }

  NumberSequance Record(boolean rec) {
    recording = rec;
    return record;
  }

  boolean IsRecursive() {
    return recursive != null;
  }

  boolean Finalize() {
    if (a >= b) {
      ResetState();
      return true;
    }
    return false;
  }

  boolean FinalizeInvertion() {
    if (b <= a) {
      ResetState();
      return true;
    }
    return false;
  }

  boolean Finished() {
    return (a >= b);
  }

  boolean FinishedInvertion() {
    return (b <= a);
  }

  void CopyResultsToVariables() {
    NumberSequance ns = this.results.Clone();
    this.variables = ns;
  }

  void CopyChangesToResults() {
    NumberSequance ns = this.changes.Clone();
    this.results = ns;
  }

  void CopyRecordToVariables() {
    if (this.record != null) {
      NumberSequance ns = this.record.Clone();
      this.variables = ns;
    }
  }

  void CopyRecordToResults() {
    if (this.record != null) {
      NumberSequance ns = this.record.Clone();
      this.results = ns;
    }
  }

  void CopyChangesToVariables() {
    NumberSequance ns = this.changes.Clone();
    this.variables = ns;
  }

  double Calc(double t) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        OneStateFunction operate = operations.First();
        Double d = results.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }
          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {

            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }
          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Calc(double t, OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a < b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }
          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Calc(double t, TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a < b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Calc() {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Calc(OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Calc(TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Sum() {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Sum(OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Sum(TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Loop() {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Loop(OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          Double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Loop(TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double LoopSum() {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double LoopSum(OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          Double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double LoopSum(TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += start;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;

          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Sum(double t) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);

          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(new Double(x));
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Loop(double t) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);

          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(new Double(x));
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double LoopSum(double t) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);

          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Sum(double t, OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Loop(double t, OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double LoopSum(double t, OneStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Sum(double t, TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append((x));
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Loop(double t, TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x = operate.f(x);
          } else {
            x = Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }
          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double LoopSum(double t, TwoStateFunction output) {
    double y = 0;
    if (operations.Length() > 0) {
      a += t;
      if (a <= b) {
        //x+=t;
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction operate = operations.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x = u.data;
          x = intergrate.f(a + x);
          if (recursive == null) {
            x += operate.f(x);
          } else {
            x += Recure(operate, x);
          }

          double r =resultant.f( x + c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((operate.f(a + u.data) - operate.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((Recure(operate, a + u.data) - Recure(operate, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          operate = operations.Next();
          if (operate == null) {
            operate = operations.First();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Inv(double t) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      b -= t;
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = (inverse.f(b + x));
          } else {
            x = InvRecure(inverse, b + x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
    }
    return y;
  }

  double Inv() {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(b + x);
          } else {
            x = InvRecure(inverse, b + x);
          }

          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          //x-=t;
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double Inv(OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(b + x);
          } else {
            x = InvRecure(inverse, b + x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double Inv(double t, OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(b + x);
          } else {
            x = InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          //x-=t;
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double Inv(TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(b + x);
          } else {
            x = InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          //x-=t;
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double Inv(double t, TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(b + x);
          } else {
            x = InvRecure(inverse, b + x);
          }

          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvSum() {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= (inverse.f(b + x));
          } else {
            x -= InvRecure(inverse, b + x);
          }

          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvSum(double t) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(b + x);
          } else {
            x -= InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          //x-=t;
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvSum(OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(b + x);
          } else {
            x -= InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvSum(double t, OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(b + x);
          } else {
            x -= InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvSum(TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(b + x);
          } else {
            x -= InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          //x-=t;
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvSum(double t, TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(b + x);
          } else {
            x -= InvRecure(inverse, b + x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          //x-=t;
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvLoopSum() {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= (inverse.f(x));
          } else {
            x -= InvRecure(inverse, x);
          }

          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          };
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvLoopSum(double t) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(x);
          } else {
            x -= InvRecure(inverse, x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvLoopSum(OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        OneStateFunction inverse = inverseOperations.Last();

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x += inverse.f(x);
          } else {
            x += InvRecure(inverse, x);
          }

          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvLoopSum(double t, OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(x);
          } else {
            x -= InvRecure(inverse, x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvLoopSum(TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x += inverse.f(x);
          } else {
            x += InvRecure(inverse, x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          //x-=t;
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvLoopSum(double t, TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {

        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x -= inverse.f(x);
          } else {
            x -= InvRecure(inverse, x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }
        }
      }
      b -= t;
    }
    return y;
  }

  double InvLoop() {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(x);
          } else {
            x = InvRecure(inverse, x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          //x-=t;
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvLoop(double t) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(x);
          } else {
            x = InvRecure(inverse, x);
          }
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          //x-=t;
          y += x;
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvLoop(OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(x);
          } else {
            x = InvRecure(inverse, x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvLoop(double t, OneStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(x);
          } else {
            x = InvRecure(inverse, x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y += output.f(x);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }

  double InvLoop(TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(x);
          } else {
            x = InvRecure(inverse, x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - start)) / start);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - start)) / start);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= start;
    }
    return y;
  }

  double InvLoop(double t, TwoStateFunction output) {
    double y = 0;
    if (inverseOperations.Length() > 0) {
      if (b >= a) {
        NumberSequance result = new NumberSequance();
        NumberSequance change = new NumberSequance();
        Double d = results.First();
        OneStateFunction inverse = inverseOperations.Last();

        for (Node<Double> u = variables.Start(); u.data != null; u = u.next) {
          double x =resultant.f(u.data);
          if (recursive == null) {
            x = inverse.f(x);
          } else {
            x = InvRecure(inverse, x);
          }
          //x-=t;
          double r = invIntergrate.f(x - c);
          result.Append(r);
          u.data = iterate.f(r);
          if (d != null) {
            if (recursive == null) {
              change.Append((inverse.f(a + u.data) - inverse.f(a + u.data - t)) / t);//(result-last result)/t
            } else {
              change.Append((InvRecure(inverse, a + u.data) - InvRecure(inverse, a + u.data - t)) / t);//(result-last result)/t
            }
          }
          y = output.f2(x, y);
          if (record != null && recording) {
            record.Append(x);
            record.n++;
          }

          inverse = inverseOperations.Prev();
          if (inverse == null) {
            inverse = inverseOperations.Last();
          }

          if (d != null) {
            d = results.Next();
          }
        }
        this.results = result;
        this.changes = change;
      }
      b -= t;
    }
    return y;
  }
}
