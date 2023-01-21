class PlanetScreen extends Screen {
  double longditude=45;
  PlanetScreen() {

    cam=new Vector();
    cam.x=radius*Math.cos(angle.x);
    cam.z=radius*Math.sin(angle.x);
    cam.y=radius*Math.cos(angle.y);
    lookat=new Vector(0, 0, 0);
    up=new Normal(0, 1, 0);
  }
  void mPressed() {

    mx=(mouseX);
    my=(mouseY);
  }
  void mReleased() {
    mx=(mouseX);
    my=(mouseY);
  }
  void mDragged() {
    double x=(mx-mouseX);
    double xs=Math.sqrt(Math.abs(x*Math.PI));
    double y=(my-mouseY);
    double ys=Math.sqrt(Math.abs(y*Math.PI));
    if (x>0) {
      x=Math.PI;
    } else if (x<0) {
      x=-Math.PI;
    }
    if (y>0) {
      y=Math.PI;
    } else if (y<0) {
      y=-Math.PI;
    }
    angle.x-=(double)(x)/width*xs;
    angle.y-=(double)(y)/height*ys;
    if (angle.x<-Math.PI) {
      angle.x+=Math.PI*2;
    } else
      if (angle.x>Math.PI) {
        angle.x-=Math.PI*2;
      }
    if (angle.y<-Math.PI) {
      angle.y+=Math.PI*2;
    } else
      if (angle.y>Math.PI) {
        angle.y-=Math.PI*2;
      }
    cam.x=radius*Math.cos(angle.x);
    cam.z=radius*Math.sin(angle.x);
    cam.y=radius*Math.cos(angle.y);
    mx=(mouseX);
    my=(mouseY);
  }
  void kPressed() {

    if (key =='w') {
      zoomUp=true;
    } else
      if (key=='s') {
        zoomDown=true;
      }
    if (keyCode==LEFT) {
      turnLeft=true;
    } else
      if (keyCode==RIGHT) {
        turnRight=true;
      }
    if (keyCode==UP) {
      turnUp=true;
    } else
      if (keyCode==DOWN) {
        turnDown=true;
      }
  }
  void kReleased() {

    if (key =='w') {
      zoomUp=false;
    } else
      if (key=='s') {
        zoomDown=false;
      }
    if (keyCode==LEFT) {
      turnLeft=false;
    } else
      if (keyCode==RIGHT) {
        turnRight=false;
      }
    if (keyCode==UP) {
      turnUp=false;
    } else
      if (keyCode==DOWN) {
        turnDown=false;
      }
  }
  void Keys() {
    if (zoomUp) {
      radius*=0.99;
      if (radius<1.25) {
        radius=1.25;
      } 

      cam.x=radius*Math.cos(angle.x);
      cam.z=radius*Math.sin(angle.x);
      cam.y=radius*Math.cos(angle.y);
    } else {
      if (zoomDown) {
        radius*=1.01;
        if (radius>25) {
          radius=25;
        }

        cam.x=radius*Math.cos(angle.x);
        cam.z=radius*Math.sin(angle.x);
        cam.y=radius*Math.cos(angle.y);
      }
    }
    double x=0;
    double y=0;
    if (turnLeft) {
      x-=Math.PI;
    }
    if (turnRight) {
      x+=Math.PI;
    }
    if (turnUp) {
      y-=Math.PI;
    }
    if (turnDown) {
      y+=Math.PI;
    }

    double xs=Math.sqrt(Math.abs(x*Math.PI));
    double ys=Math.sqrt(Math.abs(y*Math.PI));
    angle.x-=(double)(x)/width*xs;
    angle.y-=(double)(y)/height*ys;
    if (angle.x<-Math.PI) {
      angle.x+=Math.PI*2;
    } else
      if (angle.x>Math.PI) {
        angle.x-=Math.PI*2;
      }
    if (angle.y<-Math.PI) {
      angle.y+=Math.PI*2;
    } else
      if (angle.y>Math.PI) {
        angle.y-=Math.PI*2;
      }
    cam.x=radius*Math.cos(angle.x);
    cam.z=radius*Math.sin(angle.x);
    cam.y=radius*Math.cos(angle.y);
  }
  void Draw() {
    spin+=0.001f;
    longditude+=0.1f;
    if (spin>360) {
      spin-=360;
    }
    if (longditude>360) {
      longditude-=360;
    }
    float s=(float)((int)(geo.planet.size*100))/100;
    float l=(float)((int)(geo.planet.land*100))/100;
    float w=(float)((int)(geo.planet.water*100))/100;
    float t=(float)((int)(geo.planet.temp*27300-273)/PI2)/100-25;
    float a=(float)((int)(geo.planet.atmos*100))/100;
    float sc=(float)(100*geo.planet.size);
    lights();
    camera((float)cam.x*100, (float) cam.y*100, (float)cam.z*100, (float)lookat.x, (float)lookat.y, (float)lookat.z, (float)up.x, (float)up.y, (float)up.z);
    noStroke();
    pushMatrix();
    rotate(90, 0, 1, 0);
    fill(255);
    text("Size "+s+" Crust "+l+" Water "+w+" Atmospheres "+a+" Temperature "+t+"c", -200, -sc*1.25f);
    popMatrix();

    translate(0, 0, 0);

    scale(sc, sc, sc);
    rotate(spin, 0, 1, 0);
    pushMatrix();
    //translate((float)(Math.sin(spin*RADIAN)*(sc*0.75-1)), 0, (float)(Math.cos(spin*RADIAN)*(sc*0.75-1)));
    double longsin=(Math.sin(longditude*RADIAN));
    double longcos=(Math.cos(longditude*RADIAN));
    //println(longsin+":"+longcos);
    translate((float)(longcos*(1)), 0, (float)(longsin*(1)));
    scale((float)geo.planet.size, (float)geo.planet.size, (float)geo.planet.size);
    sphere(0.025f);

    popMatrix();

    DrawPlanet();
  }
}
