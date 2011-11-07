//*********************************************************************
//**      3D viewer with camera control and surface picking          **
//**              Jarek Rossignac, October 2010                      **   
//**                    (using PVectors)                             **   
//*********************************************************************
import processing.opengl.*;                // load OpenGL libraries and utilities
import javax.media.opengl.*; 
import javax.media.opengl.glu.*; 
import java.nio.*;
GL gl; 
GLU glu; 

// Global variables:
Boolean frontPick=true, translucent=true, showMesh=true, pickBack=false, showNormals=false, showVertices=false, showLabels=false; // display modes
pt F = P(0,0,0); pt T = P(); pt E = P(0,0,1000); vec U=V(0,1,0);  // focus  set with mouse when pressing '.', eye, and up vector
pt Q=P(0,0,0); vec I=V(1,0,0); vec J=V(0,1,0); vec K=V(0,0,1); // picked surface point Q and screen aligned vectors {I,J,K} set when picked
void initView() {Q=P(0,0,0); I=V(1,0,0); J=V(0,1,0); K=V(0,0,1); F = P(0,0,0); E = P(0,0,1000); U=V(0,1,0); } // declares the local frames
Mesh M = new Mesh();   // input mesh to be compressed 
Mesh R = new Mesh();   // output mesh decompressed

// SETUP
void setup() {
  size(900, 900, OPENGL); // size(500, 500, OPENGL);  
  setColors(); sphereDetail(8); 
  PFont font = loadFont("Courier-14.vlw"); textFont(font, 12);  // font for writing labels on screen
  glu= ((PGraphicsOpenGL) g).glu;  PGraphicsOpenGL pgl = (PGraphicsOpenGL) g;  gl = pgl.beginGL();  pgl.endGL();
  initView(); // declares the local frames for 3D GUI
  M.declareVectors(); R.declareVectors(); // alocates storage for points and vectors in both meshes M and R
  M.makeGrid(6); // makes a nice square mesh
         // M.loadMesh(); //  M.loadMeshOBJ(); // load mesh in my format or in OBJ (uses the next name in the fn array of names)
  M.updateON(); // computes O table and normals
  M.identifyBorderVertices();
  M.computeBox(); // computes minimax box for view setting and quantization
  F.set(M.Cbox); // sets focus on box center
  initCP(); // declare control points for curve
  }
  
// DRAW      
void draw() {  
  background(white);
  camera(E.x, E.y, E.z, F.x, F.y, F.z, U.x, U.y, U.z); // defines the view : eye, ctr, up
  vec L=U(A(V(E,F),-d(E,F),J)); directionalLight(255,255,255,L.x,L.y,L.z); // direction of light: behind and above the viewer
   
  // display model used for picking (back only when picking on the back)
  if(pickBack) {if(translucent) fill(yellow,80); else fill(255,0); noStroke(); M.showBackTriangles();}
  if(!pickBack) {
    if(translucent) {
      if((!keyPressed||key==' ')&&!mousePressed) {showCP(); showCPlabels();}
      fill(yellow,80); noStroke(); M.showBackTriangles(); fill(cyan,150); 
      if(M.showEdges) stroke(orange); 
      else noStroke(); M.showFrontTriangles();
      } 
    else {
       fill(cyan); 
       if(M.showEdges) stroke(orange); else noStroke(); 
       M.showFrontTriangles();
       }      
    }
      
  // pick point on surface    
  if (keyPressed&&key=='.') T.set(Pick()); // sets point T on the surface where the mouse points. The camera will turn toward's it when the '.' key is released
  if (keyPressed&&key=='c') M.pickc(Pick()); // picks closest corner in M
  if (keyPressed&&key=='s') M.picks(Pick()); // picks closest corner in M

  // pick control points 
  for(int i=0; i<10; i++) if (keyPressed&&key==char(i+48)) if(pickBack) CP[i].set(P(Pick(),-5,K)); else CP[i].set(P(Pick(),5,K));


  SetFrame(Q,I,J,K);  showFrame(Q,I,J,K,30);  // sets frame
  
  // view rotations
  if(keyPressed&&key==' ') {E=R(E,  PI*float(mouseX-pmouseX)/width,I,K,F); E=R(E,-PI*float(mouseY-pmouseY)/width,J,K,F); } // rotate E around F 
  if(keyPressed&&key=='Z') {E=P(E,-float(mouseY-pmouseY),K); }  //   Moves E forward/backward
  if(keyPressed&&key=='z') {E=P(E,-float(mouseY-pmouseY),K);U=R(U, -PI*float(mouseX-pmouseX)/width,I,J); }//   Moves E forward/backward and rotatees around (F,Y)

   // Display curves
  showCP(); showCPlabels(); 

 // display state
  fill(orange,200); M.showSOT(); // shoes triangle t(cc) shrunken
  M.showcc();  // display corner markers: seed sc (green),  current cc (red)
  if(showLabels) M.showLabels(); // writes IDs, but needs to be adjasted for size and orientation
  
  // display front of mesh if we were picking on the back
  if(pickBack) 
    if(translucent) {fill(cyan,150); if(M.showEdges) stroke(orange); else noStroke(); M.showFrontTriangles();} 
    else {fill(cyan); if(M.showEdges) stroke(orange); else noStroke(); M.showFrontTriangles();}
 
 // show border edges
  stroke(red); M.showBorder();
 
 // show vertices
   if(showVertices) M.showVertices();
   
 // show normals
   if(showNormals) M.showNormals();
 
 // M.showMarkedTriangles(); // shows rings around selected point
  M.showDistance();
  
  
  } // end draw
 
 // ****************** INTERRUPTS ************************* 
void mousePressed() {
  if (keyPressed&&(key=='v'||key=='d'||key=='V'||key=='D')) M.pickcOfClosestVertex(Pick()); // picks closest corner
  if (keyPressed&&key=='x') {M.pickc(Pick()); M.hide(); }// picks closest corner
  }
  
void mouseDragged() {
  if(keyPressed&&key=='v') {M.add(float(mouseX-pmouseX),I).add(-float(mouseY-pmouseY),J); } // move selected vertex in screen plane
  if(keyPressed&&key=='d') {M.add(float(mouseX-pmouseX),I).add(float(mouseY-pmouseY),K);}  // move selected vertex in X/Z screen plane
  if(keyPressed&&key=='V') {M.addROI(float(mouseX-pmouseX),I).addROI(-float(mouseY-pmouseY),J); } // move selected vertex in screen plane
  if(keyPressed&&key=='D') {M.addROI(float(mouseX-pmouseX),I).addROI(float(mouseY-pmouseY),K);}  // move selected vertex in X/Z screen plane
  }

void mouseReleased() {M.normals();}
  
void keyReleased() {
  if(key=='.') F.set(T);  // set camera focus
  U.set(M(J)); // reset camera up vector
  } 
  
void keyPressed() {
  
               // corner ops for demos
  if(key=='N') M.next();      
  if(key=='P') M.previous();
  if(key=='O') M.back();
  if(key=='L') M.left();
  if(key=='R') M.right();
  if(key=='S') M.swing();
  if(key=='U') M.unswing();
  
               // camera focus set 
  if(key=='[') F.set(M.g()); // to picked corner
  if(key==']') F.set(M.Cbox);  // center of minimax box
  if(key==';') {initView();   F.set(M.Cbox); } // reset the view
  
               // display modes
  if(key=='=') translucent=!translucent;
  if(key=='_') M.flatShading=!M.flatShading;
  if(key=='-') M.showEdges=!M.showEdges;
  if(key=='^') showVertices=!showVertices;
  if(key=='#') showLabels=!showLabels;
  if(key=='|') {showNormals=!showNormals; M.normals();}

               // archival
  if(key=='W') {M.saveMesh();}
  if(key=='G') {M.loadMesh(); M.updateON(); M.computeBox();  F.set(M.Cbox); M.identifyBorderVertices(); fni=(fni+1)%fniMax; println("Will read model "+fn[fni]);}
  
               // mesh edits, smoothing, refinement
  if(key=='b') {pickBack=true; println("picking on the back");}
  if(key=='f') {pickBack=false; println("picking on the front");}
  if(key=='/') M.flip(); // clip edge opposite to M.cc
  if(key=='B') {M.smoothen(); M.normals();}
  if(key=='Y') {M.refine();}
  if(key=='`') {M.clean();}
  if(key=='p') {
    M.clearMt();
    //M.setMt(M.computePath(M.sc,M.cc));
    /*M.setMt(M.computePath(M.retClosestCorner(CP[2]),M.retClosestCorner(CP[1])));
    M.cc=M.retClosestCorner(CP[1]);
    M.sc=M.retClosestCorner(CP[2]);
    println("cc: "+M.cc+" sc: "+M.sc);
    println("CP[1]: "+M.retClosestCorner(CP[1])+" CP[2]: "+M.retClosestCorner(CP[2]));
    */spanningTree();
  }

  if(key=='M'){
    M.makeGrid(6);
    M.clean();
  }
  

  if(key=='Q') exit();
  // M.writeCorner(); 
  } 

void spanningTree(){
  int numMarkers = 10;
  int[] visited = new int[numMarkers];
  for(int i=0; i<numMarkers; i++){
    if(CP[i].x==0 && CP[i].y==0 && CP[i].z==0)
      visited[i]=1;
  }

  while(countOnes(visited)<numMarkers){
    int i = 0;
    while(visited[i]==1){
      i++;
    }
    
    //println(i);
    visited[i]=1;

    //find closesest point CP[x] to current point CP[i] from the not visited list
    int minDist = Integer.MAX_VALUE;
    int[] minMt = new int[M.maxnt];                 // triangle markers for distance and other things   
    int[] tempMinMt = new int[M.maxnt];
    int cpMin = 0;
    for(int j=0; j<numMarkers;j++){
      if(visited[j]==0){
        int oldMinDist = minDist;
        tempMinMt = M.computePath(M.retClosestCorner(CP[i]),M.retClosestCorner(CP[j]));
        //println("comparing "+i+" and "+j);
        minDist = min(minDist,countOnes(tempMinMt));
        if(minDist != oldMinDist){
          cpMin = j;
          for(int k=0;k<M.maxnt;k++){
            minMt[k]=tempMinMt[k];
            //print(minMt[k]);
          }
        }
      }
    }
    for(int z=0;z<minMt.length;z++){
      //print(minMt[z]);
    }
    
    //visited[cpMin] = 1;
    M.addToMt(minMt);
  }
}

int countOnes(int[] x){
  int ret = 0;
  for(int i=0; i<x.length; i++){
    ret += x[i];
  }
  return ret;
}
