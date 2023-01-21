GeoSphere geo1;
GeoSphere geo2;
GeoSphere geo;
PlanetScreen pScreen;
//double au=1.5e11/6e6;//pixels
StarSystem sys;
Crono crono=new Crono();
ArrayList<Tessa> geolist=new ArrayList<Tessa>();
StarSystemScreen ssScreen;
Vector coord;
Screen cScreen;
//
void SetUpStarSystem() {

  coord=new Vector();
  coord.x=0;
  coord.y=0;
  coord.z=1000;
  double MAXX=1000;
  double MAXY=1000;
  double MAXZ=1000;
  arcRandomSeed((long)(MAXX+coord.x+(MAXY+coord.y)*MAXY+(MAXZ+coord.z)*MAXX*MAXY));
  int n0=1+(int)(arcRandom(3)+arcRandom(2));//6;
  int n1=n0+(int)(arcRandom(3)+arcRandom(2));//12;

  sys=new StarSystem(n1, 10, null, 4, 5);
  ssScreen=new StarSystemScreen();
}
String screen="Planet";
//String screen="StarSystem";
void setup() {
  size(1024, 768, P3D);
  background(0);
  SetUpStarSystem();
  geo1=new GeoSphere(5,sys);
  geo=geo1;
  pScreen=new PlanetScreen();
  geolist=geo.GetVerts(pScreen.cam, 6);
  textFont(createFont("Uroob", 24));

  if (screen=="StarSystem") {
    cScreen=ssScreen;
  } else {

    if (screen=="Planet") {
      cScreen=pScreen;
    }
  }
}
void loop() {
  //geo.Colorize(geo.planet.stain);
  //geolist=geo.GetVerts(cam,6);
}
void mousePressed() {
  cScreen.mPressed();
}
void mouseReleased() {
  cScreen.mReleased();
}
void mouseDragged() {

  cScreen.mDragged();
}
void keyPressed() {

  cScreen.kPressed();
}
void keyReleased() {

  cScreen.kReleased();
}

double RADIAN=Math.PI/180.0;
void draw() {
  background(0);
  cScreen.Keys();
  lights();
  cScreen.Draw(); 
  /*fill(255);
   scale(0.1,0.1,0.1);
   stroke(255);
   line(-300,0,0,300,0,0);
   line(0,-300,0,0,300,0);
   line(0,0,-300,0,0,300);*/
}
