class Crono {

  int epixel=250000;
  int counter=0;
  float lys=7.0;
  double scale=1.0/(lys*epixel);
  double year=(24*60*60*356);
  double day=(24*60*60);
  double week=day*7;
  double month=day*31;
  double time=week;
  double lightYear=(3.0e8)*year;
  float lspd=0.6;
}
