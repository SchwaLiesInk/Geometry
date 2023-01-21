class StellarBody extends StarSystem {
  String name;
  String dataString;
  float a;
  float b;
  double v;
  float pr;
  double rrad;
  double rad;
  double gx;
  double gy;
  float m;
  float d;
  float temp;
  StellarBody(StarSystem syst) {
    super(0, 0, syst, 1,1);
  }
}
class StarSystem extends Point {

  GImage trace;
  int stype;
  int sclass;
  boolean start=false;
  float rot=0;
  float sc=0.8+arcRandom(0.2);
  float m;
  float s1=1;//1+arcRandom(1)-arcRandom(1);
  float s2=1;//(s1+arcRandom(s1))*0.5;
  //color
  float b3;
  float g3;
  float r3;
  float sm;//star mass
  ArrayList<StellarBody> bodies=new ArrayList<StellarBody>();
  StarSystem(int nb, int mb, StarSystem sys, int starType, int starClass) {
    stype=starType;
    sclass=starClass;
    g3=(4+starType+arcRandom(2)-arcRandom(2))*32;
    r3=(10-starType+arcRandom(2)-arcRandom(2))*64;
    b3=(15-r3+arcRandom(g3)-arcRandom(2))*128;
    trace=new GImage(width*2, height*2);
    if(sys==null) {
      m=mb;
      sm=(m*m+arcRandom(m*m)+arcRandom(m*m)+arcRandom(m*m))*starType*starClass;

      int n2=(nb+(int)(arcRandom(sm/(m*m*m)))*stype*sclass)/20;
      float tm=0;

      float hz1=(float)((stype*sclass))*0.1f*0.5f;
      float hz2=(float)((stype*sclass))*0.2f*0.5f;
      for (int i=0; i<n2; i++) {

        StellarBody b=new StellarBody(this);  
        b.a=arcRandom(PI*2);
        b.m=(m*(float)(1+i)/n2+arcRandom(m*(float)(1+i)/n2)*(float)(1+i)/n2-arcRandom(m*(float)(1+i)/n2)*(float)(1+i)/n2);
        b.m*=b.m;
        tm+=b.m;
        b.pr=(float)((b.m*2)/(m+Math.sqrt(b.m)))*(1+i);

        float sg1=(float)SurfaceGravityMagnitude(b.m*6e24, b.pr*6.4e6);
        b.pr*=sg1*0.102;
        b.pr*=1+(float)(i)/n2;
        //
        sg1=(float)SurfaceGravityMagnitude(b.m*6e24, b.pr*6.4e6);
        //println(i+" sg1 "+sg1*0.102+" Mass "+b.m+" Radius "+b.pr);
        b.v=1;
        bodies.add(b);
      }
      int n3=(int)(arcRandom(sm/(m*m)));
      for (int i=0; i<n3; i++) {
        StellarBody b=new StellarBody(this);  
        b.a=arcRandom(PI*2);
        b.m=(m*(float)(1+i)/n3+arcRandom(m*(float)(1+i)/n3)*(float)(1+i)/n3-arcRandom(m*(float)(1+i)/n3)*(float)(1+i)/n3)*0.5;
        b.m*=b.m;

        tm+=b.m;
        b.pr=(float)((b.m*2)/(m+Math.sqrt(b.m)))*4;

        float sg1=(float)SurfaceGravityMagnitude(b.m*6e24, b.pr*6.4e6);
        b.pr*=sg1*0.102;
        b.pr*=1+(float)(i)/n3;
        //
        sg1=(float)SurfaceGravityMagnitude(b.m*6e24, b.pr*6.4e6);
        //println((i+1)+" sg1 "+sg1*0.102+" Mass "+b.m+" Radius "+b.pr);
        b.v=1;
        bodies.add(b);
      }

      for (int i=0; i<bodies.size (); i++) {
        int det=(int)(arcRandom(n2-i)-arcRandom(n2-i));

        if (det>0) {
          StellarBody b=bodies.get(det);
          bodies.remove(det);
          bodies.add(b);
        }
      }

      float si=(float)(crono.epixel/sm)/(tm);
      float dis=si*0.5f;
      for (int i=0; i<bodies.size (); i++) {

        StellarBody b=bodies.get(i);
        dis=b.d=((arcRandom(si)+arcRandom(2*sqrt(si*dis))))*(1+i*si);
        b.x=b.d*Math.sin(b.a);
        b.y=b.d*Math.cos(b.a);
        b.m*=(1+i*2)/(1+i*i*0.25f);//accrete
        b.pr*=(1+i*0.5)/(1+i*i*0.5f);//accrete
        b.temp=-18+((stype*sclass)+(pow(b.pr,2)))-dis*dis;//(273+16)/16
        if(b.temp<-18){b.temp=-18;}
        //
        float sg1=(float)SurfaceGravityMagnitude(b.m*6e24, b.pr*6.4e6);
        println((i+1)+" sg1 "+sg1*0.102+" Mass "+b.m+" Radius "+b.pr+" AU:"+b.d*2+" Space Temp:"+(b.temp*15.2)+"c");
      }
      float ltz=(-18+((stype*sclass)+(pow(hz1,2)))-hz2*hz2)*15.2;
      float utz=(-18+((stype*sclass)+(pow(hz2,2)))-hz1*hz1)*15.2;
      println("solar "+sm/3000+" effect"+si+" on "+tm+" HZ:"+hz1+"->"+hz2+" HZ Space Temp:"+utz+"/"+ltz);
      b3*=16;
      r3*=16;
      g3*=16;
    }
  }
  void Draw() {
    float xs=0;//width*0.5f;
    float ys=0;//height*0.5f;
    pushMatrix();
    scale(2.0, 2.0, 2.0);
    trace.Draw3D();
    popMatrix();

    scale(0.5, 0.5, 0.5);
    fill(255);
    //text("At "+((1.0-scale)*lys), 100, 60);

    float na=(float)(0.5/(crono.epixel*PI*8));
    start=true;
    crono.scale=1.0;
    /*if (scale<1.0) {
     scale+=(lys*epixel*time)/((lightYear)*lspd);
     trak=32;
     } else {
     trak=8;
     start=true;
     }*/
    sc=(float)(crono.scale);
    sc*=sc*sc*sc*sc;
    sc*=45;
    //double lax=0;
    //double ax=0;
    //double lay=0;
    //double ay=0;

    for (int i=0; i<bodies.size (); i++) {
      StellarBody b1=bodies.get(i);    
      double srx1=(-b1.x)*1.5e11;
      double sry1=(-b1.y)*1.5e11;
      double stx1=(double)((srx1*srx1)*3e8);
      double sty1=(double)((sry1*sry1)*3e8);
      //
      b1.gx=(uG*2*(b1.m*6e24+sm*6e24))/(stx1);
      b1.gy=(uG*2*(b1.m*6e24+sm*6e24))/(sty1);
      for (int j=0; j<bodies.size (); j++) {
        StellarBody b2=bodies.get(j);      
        if (b1!=b2) {
          double x1=((b2.x-b1.x)*1.5e11);
          double y1=((b2.y-b1.y)*1.5e11);

          double xt1=(double)((x1*x1)*3e8);
          double yt1=(double)((y1*y1)*3e8);
          b1.gx+=(uG*2*(b1.m*6e24+b2.m*6e24))/(xt1);
          b1.gy+=(uG*2*(b1.m*6e24+b2.m*6e24))/(yt1);
        }
      }
      b1.b=(float)Math.atan(b1.gx/b1.gy);
    }
    int l=0;

    beginShape(LINES);
    for (int k=0; k<64; k+=l) {
      stroke(255, 255, 255, 64-k);
      vertex(xs+sc*(1+k), ys-(1+k)*sc);
      vertex(xs+sc*(1+k), ys+(1+k)*sc);//au
      vertex(xs-sc*(1+k), ys-(1+k)*sc);
      vertex(xs-sc*(1+k), ys+(1+k)*sc);//au
      vertex(xs-sc*(1+k), ys+(1+k)*sc);
      vertex(xs+sc*(1+k), ys+(1+k)*sc);//au
      vertex(xs-sc*(1+k), ys-(1+k)*sc);
      vertex(xs+sc*(1+k), ys-(1+k)*sc);//au
      l+=1;
    }
    endShape();
    float ra=(float)(sm/(m*m)*sc*0.0005*sclass);
    int sstype=(stype*sclass);
    for (int k=0; k<16; k++) {
      float kr=k*ra;
      //float kd=kr*0.001;
      float kc=(255-k*k)*2;
      if (kc<1)kc=1;
      noStroke();
      //stroke((int)(r3*kc)%255, (int)(g3*kc)%255, (int)(b3*kc)%255,4);
      fill((int)(r3+kc)/sstype, (int)(g3+kc)/sstype, (int)(b3+kc)/sstype, kc );
      translate(xs, 0, ys);
      sphere(kr);
    }
    //fill(255,255, 255, 32 );
    //sphere(ra*32);
    //fill(255,255, 255, 16 );
    //sphere(ra*64);

    for (int i=0; i<bodies.size (); i++) {
      StellarBody b=bodies.get(i);

      float gx1=(float)(b.gx*crono.time);
      float gy1=(float)(b.gy*crono.time);
      //float ara=(float)(Math.sqrt(x*x+y*y+3e-8));
      double gr=Math.sqrt(b.gx*b.gx+b.gy*b.gy+3e-8)*crono.time;//gravocity
      float ar=(float)(Math.sqrt(b.x*b.x+b.y*b.y+3e-8));
      //double t=(Math.atan(b.gx/b.gy));
      //double s=Math.sin(t);
      //double c=Math.cos(t+PI);
      double ga=Math.sqrt(b.gx*b.gx+b.gy*b.gy+3e-8);
      double tps=(((crono.time)/((PI*2*ga)*(3e8))));
      double sv=crono.time/tps;
      b.v=gr+sv/b.v;
      float rc1=(sqrt(ar));
      float va=(float)(na*b.v*crono.time/(6e6*ar*ar*PI*8))*crono.epixel*0.5/(stype*sclass);
      b.a+=va;//;
      b.rad=b.a;
      
      b.x+=(8*PI*sv*(crono.time/3e8))/(ar*ar*1.5e11*PI*2)*Math.sin(b.a);//em;
      b.y+=(8*PI*sv*(crono.time/3e8))/(ar*ar*1.5e11*PI*2)*Math.cos(b.a);//em;
      b.x+=(8*PI*sv*(crono.time/3e8))/(ar*1.5e11*PI*2)*Math.sin(b.a);//-b.gx*time;
      b.y+=(8*PI*sv*(crono.time/3e8))/(ar*1.5e11*PI*2)*Math.cos(b.a);//-b.gy*time;

      b.rrad=Math.sqrt(b.x*b.x+b.y*b.y+3e-8)*(stype*sclass);

      //float rot2=(float)(b.b);
      float yn1=(float)((b.rrad)*Math.sin(b.a));
      float yn2=(float)((b.rrad)*Math.cos(b.a));
      float grav1=(float)Math.sqrt(SurfaceGravityMagnitude(b.m*6e24, b.pr*6.4e6))*0.25;
      float g1=(rc1*2)*128;
      float r1=(rc1)*128;
      float b1=(2.0/rc1)*128;
      if (r1>255) {
        r1=255;
      }
      if (g1>255) {
        g1=255;
      }
      if (b1>255) {
        b1=255;
      }
      float sc1=(float)(sc*s1);
      float sc2=(float)(sc*s2);

      stroke(255-(int)(r1*rc1)%255, 255-(int)(g1*rc1)%255, 255-(int)(b1*rc1)%255, 255);

      point(xs+(yn1+gx1)*sc1, ys+(yn2+gy1)*sc2);
      noFill();
      rect(xs+(yn1+gx1)*sc1-grav1*sc*0.5, ys+(yn2+gy1)*sc2-grav1*sc*0.5, grav1*sc, grav1*sc);
      fill(255-(int)(r1*rc1)%255, 255-(int)(g1*rc1)%255, 255-(int)(b1*rc1)%255, 32);
      noStroke();
      pushMatrix();
      translate(xs+(yn1+gx1)*sc1, ys+(yn2+gy1)*sc2, 0);
      sphere(b.m*0.0025*sc);
      popMatrix();
      fill(255-(int)(r1*rc1)%255, 255-(int)(g1*rc1)%255, 255-(int)(b1*rc1)%255, 64);
      pushMatrix();
      translate(xs+(yn1+gx1)*sc1, ys+(yn2+gy1)*sc2, 0);
      sphere((float)(gr*sc*0.0025));
      popMatrix();

      if (start) {

        trace.AddPixel(trace.bits.width*0.5+(yn1+gx1)*sc1*0.25, trace.bits.height*0.5+(yn2+gy1)*sc2*0.25, 128, 255-(int)(r1*rc1)%255, 255-(int)(g1*rc1)%255, 255-(int)(b1*rc1)%255);
      }
    }
    crono.counter++;
    if (start && crono.counter>2) {
      trace.ClearAdd(new Pigment(16, 32, 64, 255), 1500, 1.0f);
      crono.counter=0;
    }
    trace.Update();
  }
}
