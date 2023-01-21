abstract class Screen {

  Vector cam=new Vector(0, 0, 10);
  Vector angle=new Vector(0, Math.PI*0.5, 0);
  Vector lookat=new Vector(0, 0, 0);
  Normal up=new Normal(0, 1, 0);

  double radius=3;
  float spin=0;
  int mx;  
  int my;  
  boolean zoomUp=false;
  boolean zoomDown=false; 
  boolean turnLeft=false;
  boolean turnRight=false; 
  boolean turnUp=false;
  boolean turnDown=false; 
  abstract void mPressed();
  abstract void mReleased();
  abstract void mDragged();
  abstract void kPressed();
  abstract void kReleased();
  abstract void Keys();
  abstract void Draw();
}
