// *** MARKERS ON THE MESH ****
pt [] CP = new pt [10];
void initCP() {
  for (int i=0; i<CP.length; i++) CP[i]=P(); 
  //float r=100; for (int i=1; i<7; i++) {CP[i].x=r*cos(TWO_PI*i/6); CP[i].y=r*sin(TWO_PI*i/6);}
  float r=100; for (int i=1; i<7; i++) {CP[i].x=0; CP[i].y=0;}
  } // creates points
void showCP() {noStroke(); for (int i=0; i<CP.length; i++) {fill(ramp(i,CP.length)); show(CP[i],10); }} // shows control points
void showCPlabels() {noStroke(); fill(black); for (int i=0; i<CP.length; i++) {fill(ramp(i,CP.length)); show(CP[i],str(i),V(10,I,10,J,2,K));} } // shows IDs
   
  
