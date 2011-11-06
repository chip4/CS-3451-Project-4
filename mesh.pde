// CORNER TABLE FOR TRIANGLE MESHES by Jarek Rosignac
// Last edited October, 2011
// example meshesshowShrunkOffsetT
String [] fn= {"bunny.vts","horse.vts","torus.vts","flat.vts","tet.vts","fandisk.vts","squirrel.vts","venus.vts"};
int fni=0; int fniMax=fn.length; // file names for loading meshes

//========================== class MESH ===============================
class Mesh {
//  ==================================== Internal variables ====================================
 // max sizes, counts, selected corners
 int maxnv = 45000;                         //  max number of vertices
 int maxnt = maxnv*2;                       // max number of triangles
 int nv = 0;                              // current  number of vertices
 int nt = 0;                   // current number of triangles
 int nc = 0;                                // current number of corners (3 per triangle)
 int cc=0, pc=0, sc=0;                      // current, previous, saved corners
 
 // primary tables
 int[] V = new int [3*maxnt];               // V table (triangle/vertex indices)
 int[] O = new int [3*maxnt];               // O table (opposite corner indices)
 pt[] G = new pt [maxnv];                   // geometry table (vertices)


vec[] Nv = new vec [maxnv];                 // vertex normals or laplace vectors
vec[] Nt = new vec [maxnt];                // triangles normals

 
 // auxiliary tables for bookkeeping
 int[] cm = new int[3*maxnt];               // corner markers: 
 int[] vm = new int[3*maxnt];               // vertex markers: 0=not marked, 1=interior, 2=border, 3=non manifold
 int[] tm = new int[3*maxnt];               // triangle markers: 0=not marked, 

 // other tables
 int[] Mv = new int[maxnv];                  // vertex markers
 int [] Valence = new int [maxnv];          // vertex valence (count of incident triangles)

 int[] Mt = new int[maxnt];                 // triangle markers for distance and other things   
 boolean [] VisitedT = new boolean [maxnt];  // triangle visited
 boolean[] visible = new boolean[maxnt];    // set if triangle visible

 int[] W = new int [3*maxnt];               // mid-edge vertex indices for subdivision (associated with corner opposite to edge)

 pt[] G2 = new pt [maxnv]; //2008-03-06 JJ misc
 boolean [] Border = new boolean [maxnv];   // vertex is border
 boolean [] VisitedV = new boolean [maxnv];  // vertex visited
 int r=2;                                // radius of spheres for displaying vertices
 float [] distance = new float [maxnv];
 // geodesic
 boolean  showPath=false, showDistance=false;  
 boolean[] P = new boolean [3*maxnt];       // marker of corners in a path to parent triangle
 int[] Distance = new int[maxnt];           // triangle markers for distance fields 
 int[] SMt = new int[maxnt];                // sum of triangle markers for isolation
 int prevc = 0;                             // previously selected corner
 int rings=2;                           // number of rings for colorcoding


 // box
 pt Cbox = new pt(width/2,height/2,0);                   // mini-max box center
 float rbox=1000;                                        // half-diagonal of enclosing box

 // rendering modes
 Boolean flatShading=false, showEdges=true;  // showEdges shoes edges as gaps. Smooth shading works only when !showEdge

//  ==================================== INIT, CREATE, COPY ====================================
 Mesh() {}

 void declareVectors() {
   for (int i=0; i<maxnv; i++) {G[i]=P(); Nv[i]=V();};   // init vertices and normals
   for (int i=0; i<maxnt; i++) Nt[i]=V();       // init triangle normals and skeleton lab els
   }

 void resetCounters() {nv=0; nt=0; nc=0;}

 void makeGrid (int w) { // make a 2D grid of w x w vertices
   for (int i=0; i<w; i++) {for (int j=0; j<w; j++) { G[w*i+j].set(height*.8*j/(w-1)+height/10,height*.8*i/(w-1)+height/10,0);}}    
   for (int i=0; i<w-1; i++) {for (int j=0; j<w-1; j++) {                  // define the triangles for the grid
     V[(i*(w-1)+j)*6]=i*w+j;       V[(i*(w-1)+j)*6+2]=(i+1)*w+j;       V[(i*(w-1)+j)*6+1]=(i+1)*w+j+1;
     V[(i*(w-1)+j)*6+3]=i*w+j;     V[(i*(w-1)+j)*6+5]=(i+1)*w+j+1;     V[(i*(w-1)+j)*6+4]=i*w+j+1;}; };
   nv = w*w;
   nt = 2*(w-1)*(w-1); 
   nc=3*nt;  }

 void resetMarkers() { // reset the seed and current corner and the markers for corners, triangles, and vertices
   cc=0; pc=0; sc=0;
   for (int i=0; i<nv; i++) vm[i]=0;
   for (int i=0; i<nc; i++) cm[i]=0;
   for (int i=0; i<nt; i++) tm[i]=0;
   for (int i=0; i<nt; i++) visible[i]=true;
   }
 
 int addVertex(pt P) { G[nv].set(P); nv++; return nv-1;};
 int addVertex(float x, float y, float z) { G[nv].x=x; G[nv].y=y; G[nv].z=z; nv++; return nv-1;};
  
 void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; visible[nt++]=true;} // adds a triangle

 void updateON() {computeO(); normals(); } // recomputes O and normals

  // ============================================= CORNER OPERATORS =======================================
    // operations on a corner
  int t (int c) {return int(c/3);};              // triangle of corner    
  int n (int c) {return 3*t(c)+(c+1)%3;};        // next corner in the same t(c)    
  int p (int c) {return n(n(c));};               // previous corner in the same t(c)  
  int v (int c) {return V[c] ;};                 // id of the vertex of c             
  int o (int c) {return O[c];};                  // opposite (or self if it has no opposite)
  int l (int c) {return o(n(c));};               // left neighbor (or next if n(c) has no opposite)                      
  int r (int c) {return o(p(c));};               // right neighbor (or previous if p(c) has no opposite)                    
  int s (int c) {return n(l(c));};               // swings around v(c) or around a border loop
  int u (int c) {return p(r(c));};               // unswings around v(c) or around a border loop
  int c (int t) {return t*3;}                    // first corner of triangle t
  boolean b (int c) {return O[c]==c;};           // if faces a border (has no opposite)
  boolean vis(int c) {return visible[t(c)]; };   // true if tiangle of c is visible

    // operations on the selected corner cc
  int t() {return t(cc); }
  int n() {return n(cc); }        Mesh next() {pc=cc; cc=n(cc); return this;};
  int p() {return p(cc); }        Mesh previous() {pc=cc; cc=p(cc); return this;};
  int v() {return v(cc); }
  int o() {return o(cc); }        Mesh back() {if(!b(cc)) {pc=cc; cc=o(cc);}; return this;};
  boolean b() {return b(cc);}             
  int l() {return l(cc);}         Mesh left() {next(); back(); return this;}; 
  int r() {return r(cc);}         Mesh right() {previous(); back(); return this;};
  int s() {return s(cc);}         Mesh swing() {left(); next();  return this;};
  int u() {return u(cc);}         Mesh unswing() {right(); previous();  return this;};

    // geometry for corner c
  pt g (int c) {return G[v(c)];};                // shortcut to get the point of the vertex v(c) of corner c
  pt cg(int c) {pt cPt = P(g(c),.3,triCenter(t(c)));  return(cPt); };   // computes point at corner
  pt corner(int c) {return P(g(c),triCenter(t(c)));   };   // returns corner point

    // normals fot t(c) (must be precomputed)
  vec Nv (int c) {return(Nv[V[c]]);}; vec Nv() {return Nv(cc);}            // shortcut to get the normal of v(c) 
  vec Nt (int c) {return(Nt[t(c)]);}; vec Nt() {return Nt(cc);}            // shortcut to get the normal of t(c) 

    // geometry for corner cc
  pt g() {return g(cc);}            // shortcut to get the point of the vertex v(c) of corner c
  void setG(pt P) {G[v(cc)].set(P);} // moves vertex of c to P

     // debugging prints
  void writeCorner (int c) {println("cc="+cc+", n="+n(cc)+", p="+p(cc)+", o="+o(cc)+", v="+v(cc)+", t="+t(cc)+"."+", nt="+nt+", nv="+nv ); }; 
  void writeCorner () {writeCorner (cc);}
  void writeCorners () {for (int c=0; c<nc; c++) {println("T["+c+"]="+t(c)+", visible="+visible[t(c)]+", v="+v(c)+",  o="+o(c));};}

// ============================================= MESH MANIPULATION =======================================
  // pick corner closest to point X
  void pickcOfClosestVertex (pt X) {for (int b=0; b<nc; b++) if(d(X,g(b))<d(X,g(cc))) {cc=b; pc=b; } } // picks corner of closest vertex to X
  void pickc (pt X) {for (int b=0; b<nc; b++) if(d(X,cg(b))<d(X,cg(cc))) {cc=b; pc=b; } } // picks closest corner to X
  void picksOfClosestVertex (pt X) {for (int b=0; b<nc; b++) if(d(X,g(b))<d(X,g(sc))) {sc=b;} } // picks corner of closest vertex to X
  void picks (pt X) {for (int b=0; b<nc; b++) if(d(X,cg(b))<d(X,cg(sc))) {sc=b;} } // picks closest corner to X

  boolean ptInTriangle(pt P, int t){
    //println("in ptInTriangle");
    pt A = G[v(c(t))];
    pt B = G[v(n(c(t)))];
    pt C = G[v(n(n(c(t))))];

    float areaT = areaOfT(A,B,C);
    float epsilon = 200;

    float areaTotal = areaOfT(P,B,C) + areaOfT(A,P,C) + areaOfT(A,B,P);
    //println("areaT: " + areaT + " areaTotal: "+areaTotal);
    //println(areaT-areaTotal);

    //println((areaT-epsilon <= areaTotal)&&(areaT+epsilon >= areaTotal));
    return (areaT-epsilon <= areaTotal)&&(areaT+epsilon >= areaTotal);

  } 
  /*float sumAreas(pt P, int t){
    pt A = G[v(c(t))];
    pt B = G[v(n(c(t)))];
    pt C = G[v(n(n(c(t))))];

    float areaT = areaOfT(A,B,C);
    float epsilon = areaT/5;

    float areaTotal = areaOfT(P,B,C) + areaOfT(A,P,C) + areaOfT(A,B,P);

    return abs(areaT-areaTotal);
  }*/
  float areaOfT(pt p1, pt p2, pt p3){
    return .5*sqrt(
        sq(det(p1.x,p2.x,p3.x,p1.y,p2.y,p3.y)) +
        sq(det(p1.y,p2.y,p3.y,p1.z,p2.z,p3.z)) +
        sq(det(p1.z,p2.z,p3.z,p1.x,p2.x,p3.x)) 
      );
  }
  float det(float p1, float p2, float p3, float p4, float p5, float p6){
    return (p1*(p5-p6)-p2*(p4-p6)+p3*(p4-p5));
  }
  int retClosestCorner (pt X) {
    int ret = 0;
    for (int b=0; b<nc; b++) 
      if(d(X,g(b))<d(X,g(ret))) {
        //if(sumAreas(X,t(b)) < sumAreas(X,t(ret))){
        if(ptInTriangle(X,t(b))){
          ret=b; 
        }
        //pc=b; 
      } 
    return ret;
  } // returns corner of closest vertex to X

  // move the vertex of a corner
  void setG(int c, pt P) {G[v(c)].set(P);}       // moves vertex of c to P
  Mesh add(int c, vec V) {G[v(c)].add(V); return this;}             // moves vertex of c to P
  Mesh add(int c, float s, vec V) {G[v(c)].add(s,V); return this;}   // moves vertex of c to P
  Mesh add(vec V) {G[v(cc)].add(V); return this;} // moves vertex of c to P
  Mesh add(float s, vec V) {G[v(cc)].add(s,V); return this;} // moves vertex of c to P
  void move(int c) {g(c).add(pmouseY-mouseY,Nv(c));}
  void move(int c, float d) {g(c).add(d,Nv(c));}
  void move() {move(cc); normals();}

  Mesh addROI(float s, vec V) { return addROI(64,s,V);}
  Mesh addROI(int d, float s, vec V) {
     float md=setROI(d); 
     for (int c=0; c<nc; c++) if(!VisitedV[v(c)]&&(Mv[v(c)]!=0))  G[v(c)].add(s*(1.-distance[v(c)]/md),V);   // moves ROI
     smoothROI();
     setROI(d*2); // marks ROI of d rings
     smoothROI(); smoothROI();
     return this;
     }   

  void tuckROI(float s) {for (int i=0; i<nv; i++) if (Mv[i]!=0) G[i].add(s,Nv[i]); };  // displaces each vertex by a fraction s of its normal
  void smoothROI() {computeLaplaceVectors(); tuckROI(0.5); computeLaplaceVectors(); tuckROI(-0.5);};
  
float setROI(int n) { // marks vertices and triangles at a graph distance of maxr
  float md=0;
  int tc=0; // triangle counter
  int r=1; // ring counter
  for(int i=0; i<nt; i++) {Mt[i]=0;};  // unmark all triangles
  Mt[t(cc)]=1; tc++;                   // mark t(cc)
  for(int i=0; i<nv; i++) {Mv[i]=0;};  // unmark all vertices
  while ((tc<nt)&&(tc<n)) {  // while not finished
     for(int i=0; i<nc; i++) {if ((Mv[v(i)]==0)&&(Mt[t(i)]==r)) {Mv[v(i)]=r; distance[v(i)]=d(g(cc),g(i)); md = max(md,distance[v(i)]); };};  // mark vertices of last marked triangles
     for(int i=0; i<nc; i++) {if ((Mt[t(i)]==0)&&(Mv[v(i)]==r)) {Mt[t(i)]=r+1; tc++;};}; // mark triangles incident on last marked vertices
     r++; // increment ring counter
     };
  rings=r;
  return md;
  }

 //  ==========================================================  HIDE TRIANGLES ===========================================
void markRings(int maxr) { // marks vertices and triangles at a graph distance of maxr
  int tc=0; // triangle counter
  int r=1; // ring counter
  for(int i=0; i<nt; i++) {Mt[i]=0;};  // unmark all triangles
  Mt[t(cc)]=1; tc++;                   // mark t(cc)
  for(int i=0; i<nv; i++) {Mv[i]=0;};  // unmark all vertices
  while ((tc<nt)&&(r<=maxr)) {  // while not finished
     for(int i=0; i<nc; i++) {if ((Mv[v(i)]==0)&&(Mt[t(i)]==r)) {Mv[v(i)]=r;};};  // mark vertices of last marked triangles
     for(int i=0; i<nc; i++) {if ((Mt[t(i)]==0)&&(Mv[v(i)]==r)) {Mt[t(i)]=r+1; tc++;};}; // mark triangles incident on last marked vertices
     r++; // increment ring counter
     };
  rings=r; // sets ring variable for rendring?
  }

void hide() {visible[t(cc)]=false;}

// ============================================= GEOMETRY =======================================
  // enclosing box
  void computeBox() { // computes center Cbox and half-diagonal Rbox of minimax box
    pt Lbox =  P(G[0]);  pt Hbox =  P(G[0]);
    for (int i=1; i<nv; i++) { 
      Lbox.x=min(Lbox.x,G[i].x); Lbox.y=min(Lbox.y,G[i].y); Lbox.z=min(Lbox.z,G[i].z);
      Hbox.x=max(Hbox.x,G[i].x); Hbox.y=max(Hbox.y,G[i].y); Hbox.z=max(Hbox.z,G[i].z); 
      };
    Cbox.set(P(Lbox,Hbox));  rbox=d(Cbox,Hbox); 
    };

// ============================================= O TABLE CONSTRUCTION =========================================
  void computeOnaive() {                        // sets the O table from the V table, assumes consistent orientation of triangles
    resetCounters();
    for (int i=0; i<3*nt; i++) {O[i]=i;};  // init O table to -1: has no opposite (i.e. is a border corner)
    for (int i=0; i<nc; i++) {  for (int j=i+1; j<nc; j++) {       // for each corner i, for each other corner j
        if( (v(n(i))==v(p(j))) && (v(p(i))==v(n(j))) ) {O[i]=j; O[j]=i;};};};}// make i and j opposite if they match         

  void computeO() {
    resetMarkers(); 
    int val[] = new int [nv]; for (int v=0; v<nv; v++) val[v]=0;  for (int c=0; c<nc; c++) val[v(c)]++;   //  valences
    int fic[] = new int [nv]; int rfic=0; for (int v=0; v<nv; v++) {fic[v]=rfic; rfic+=val[v];};  // head of list of incident corners
    for (int v=0; v<nv; v++) val[v]=0;   // valences wil be reused to track how many incident corners were encountered for each vertex
    int [] C = new int [nc]; for (int c=0; c<nc; c++) C[fic[v(c)]+val[v(c)]++]=c;  // vor each vertex: the list of val[v] incident corners starts at C[fic[v]]
    for (int c=0; c<nc; c++) O[c]=c;    // init O table to -1 meaning that a corner has no opposite (i.e. faces a border)
    for (int v=0; v<nv; v++)             // for each vertex...
       for (int a=fic[v]; a<fic[v]+val[v]-1; a++) for (int b=a+1; b<fic[v]+val[v]; b++)  { // for each pair (C[a],C[b[]) of its incident corners
          if (v(n(C[a]))==v(p(C[b]))) {O[p(C[a])]=n(C[b]); O[n(C[b])]=p(C[a]); }; // if C[a] follows C[b] around v, then p(C[a]) and n(C[b]) are opposite
          if (v(n(C[b]))==v(p(C[a]))) {O[p(C[b])]=n(C[a]); O[n(C[a])]=p(C[b]); };        
        };                
     }

// ============================================= DISPLAY CORNERS and LABELS =============================
  void showCorner(int c, float r) {show(cg(c),r); };   // renders corner c as small ball
  
  void showcc(){noStroke(); fill(green); showCorner(sc,4); /* fill(green); showCorner(pc,5); */ fill(red); showCorner(cc,6); } // displays corner markers
  
  void showLabels() { // displays IDs of corners, vertices, and triangles
   fill(black); 
   for (int i=0; i<nv; i++) {show(G[i],"v"+str(i),V(10,Nv[i])); }; 
   for (int i=0; i<nc; i++) {show(corner(i),"c"+str(i),V(10,Nt[i])); }; 
   for (int i=0; i<nt; i++) {show(triCenter(i),"t"+str(i),V(10,Nt[i])); }; 
   noFill();
   }

// ============================================= DISPLAY VERTICES =======================================
  void showVertices () {
    noStroke(); noSmooth(); 
    for (int v=0; v<nv; v++)  {
      if (vm[v]==0) fill(yellow,150);
      if (vm[v]==1) fill(red,150);
      if (vm[v]==2) fill(green,150);
      if (vm[v]==3) fill(blue,150);
      if(Border[v]) fill(magenta,150);
       show(G[v],r);  
      }
    noFill();
    }

// ============================================= DISPLAY EDGES =======================================
  void showBorder() {for (int c=0; c<nc; c++) {if (b(c)) {drawEdge(c);}; }; };         // draws all border edges
  void showEdges () {for(int c=0; c<nc; c++) drawEdge(c); };  
  void drawEdge(int c) {show(g(p(c)),g(n(c))); };  // draws edge of t(c) opposite to corner c

// ============================================= DISPLAY TRIANGLES =======================================
  // displays triangle if marked as visible using flat or smooth shading (depending on flatShading variable
  void shade(int t) { // displays tiangle t if visible
    if(visible[t])  
      if(flatShading) {beginShape(); vertex(g(3*t)); vertex(g(3*t+1)); vertex(g(3*t+2));  endShape(CLOSE); }
      else {beginShape(); normal(Nv[v(3*t)]); vertex(g(3*t)); normal(Nv[v(3*t+1)]); vertex(g(3*t+1)); normal(Nv[v(3*t+2)]); vertex(g(3*t+2));  endShape(CLOSE); }; 
    }
  
  // display shrunken and offset triangles
  void showShrunkT(int t, float e) {if(visible[t]) showShrunk(g(3*t),g(3*t+1),g(3*t+2),e);}
  void showSOT(int t) {if(visible[t]) showShrunkOffsetT(t,1,1);}
  void showSOT() {if(visible[t(cc)]) showShrunkOffsetT(t(cc),1,1);}
  void showShrunkOffsetT(int t, float e, float h) {if(visible[t]) showShrunkOffset(g(3*t),g(3*t+1),g(3*t+2),e,h);}
  void showShrunkT() {int t=t(cc); if(visible[t]) showShrunk(g(3*t),g(3*t+1),g(3*t+2),2);}
  void showShrunkOffsetT(float h) {int t=t(cc); if(visible[t]) showShrunkOffset(g(3*t),g(3*t+1),g(3*t+2),2,h);}

  // display front and back triangles shrunken if showEdges  
  Boolean frontFacing(int t) {return !cw(E,g(3*t),g(3*t+1),g(3*t+2)); } 
  void showFrontTriangles() {for(int t=0; t<nt; t++) if(frontFacing(t)) {if(showEdges) showShrunkT(t,1); else shade(t);}};  
  void showBackTriangles() {for(int t=0; t<nt; t++) if(!frontFacing(t)) if(showEdges) showShrunkT(t,1); else shade(t);};  
  void showAllTriangles() {for(int t=0; t<nt; t++) if(showEdges) showShrunkT(t,1); else shade(t);};  
  void showMarkedTriangles() {for(int t=0; t<nt; t++) if(visible[t]&&Mt[t]!=0) {fill(ramp(Mt[t],rings)); showShrunkOffsetT(t,1,1); }};  

//  ==========================================================  PROCESS EDGES ===========================================
  // FLIP 
  void flip(int c) {      // flip edge opposite to corner c, FIX border cases
    if (b(c)) return;
      V[n(o(c))]=v(c); V[n(c)]=v(o(c));
      int co=o(c); 
      O[co]=r(c); 
      if(!b(p(c))) O[r(c)]=co; 
      if(!b(p(co))) O[c]=r(co); 
      if(!b(p(co))) O[r(co)]=c; 
      O[p(c)]=p(co); O[p(co)]=p(c);  
    }
  void flip() {flip(cc); pc=cc; cc=p(cc);}

  void flipWhenLonger() {for (int c=0; c<nc; c++) if (d(g(n(c)),g(p(c)))>d(g(c),g(o(c)))) flip(c); } 

  int cornerOfShortestEdge() {  // assumes manifold
    float md=d(g(p(0)),g(n(0))); int ma=0;
    for (int a=1; a<nc; a++) if (vis(a)&&(d(g(p(a)),g(n(a)))<md)) {ma=a; md=d(g(p(a)),g(n(a)));}; 
    return ma;
    } 
  void findShortestEdge() {cc=cornerOfShortestEdge();  } 

//  ========================================================== PROCESS  TRIANGLES ===========================================
 pt triCenter(int i) {return P( G[V[3*i]], G[V[3*i+1]], G[V[3*i+2]] ); };  
 pt triCenter() {return triCenter(t());}  // computes center of triangle t(i) 
 void writeTri (int i) {println("T"+i+": V = ("+V[3*i]+":"+v(o(3*i))+","+V[3*i+1]+":"+v(o(3*i+1))+","+V[3*i+2]+":"+v(o(3*i+2))+")"); };

 
   
//  ==========================================================  NORMALS ===========================================
void normals() {computeTriNormals(); computeVertexNormals(); }
void computeValenceAndResetNormals() {      // caches valence of each vertex
  for (int i=0; i<nv; i++) {Nv[i]=V();  Valence[i]=0;};  // resets the valences to 0
  for (int i=0; i<nc; i++) {Valence[v(i)]++; };
  }
vec triNormal(int t) { return N(V(g(3*t),g(3*t+1)),V(g(3*t),g(3*t+2))); };  
void computeTriNormals() {for (int i=0; i<nt; i++) {Nt[i].set(triNormal(i)); }; };             // caches normals of all tirangles
void computeVertexNormals() {  // computes the vertex normals as sums of the normal vectors of incident tirangles scaled by area/2
  for (int i=0; i<nv; i++) {Nv[i].set(0,0,0);};  // resets the valences to 0
  for (int i=0; i<nc; i++) {Nv[v(i)].add(Nt[t(i)]);};
  for (int i=0; i<nv; i++) {Nv[i].normalize();};            };
void showVertexNormals() {for (int i=0; i<nv; i++) show(G[i],V(10*r,Nv[i]));  };
void showTriNormals() {for (int i=0; i<nt; i++) show(triCenter(i),V(10*r,U(Nt[i])));  };
void showNormals() {if(flatShading) showTriNormals(); else showVertexNormals(); }

// ============================================================= SMOOTHING ============================================================
void identifyBorderVertices() {
  for (int v=0; v<nv; v++) Border[v]=false;
  for (int c=0; c<nc; c++) if(b(c)) Border[v(n(c))]=true;
  }
void computeLaplaceVectors() {  // computes the vertex normals as sums of the normal vectors of incident triangles scaled by triangle area*2
  computeValenceAndResetNormals();
  for (int i=0; i<3*nt; i++) {Nv[v(p(i))].add(V(g(p(i)),g(n(i))));};
  for (int i=0; i<nv; i++) {Nv[i].div(Valence[i]);};                         };
void tuck(float s) {for (int i=0; i<nv; i++) if(!Border[i]) G[i].add(s,Nv[i]); };  // displaces each vertex by a fraction s of its normal
void smoothen() {normals(); computeLaplaceVectors(); tuck(0.6); computeLaplaceVectors(); tuck(-0.6);};

// ============================================================= SUBDIVISION ============================================================
int w (int c) {return(W[c]);};               // temporary indices to mid-edge vertices associated with corners during subdivision

void splitEdges() {            // creates a new vertex for each edge and stores its ID in the W of the corner (and of its opposite if any)
  for (int i=0; i<3*nt; i++) {  // for each corner i
    if(b(i)) {G[nv]=P(g(n(i)),g(p(i))); W[i]=nv++;}
    else {if(i<o(i)) {G[nv]=P(g(n(i)),g(p(i))); W[o(i)]=nv; W[i]=nv++; }; }; }; } // if this corner is the first to see the edge
  
void bulge() {              // tweaks the new mid-edge vertices according to the Butterfly mask
  for (int i=0; i<3*nt; i++) {
    if((!b(i))&&(i<o(i))) {    // no tweak for mid-vertices of border edges
     if (!b(p(i))&&!b(n(i))&&!b(p(o(i)))&&!b(n(o(i))))
      {G[W[i]].add(0.25,V(P(P(g(l(i)),g(r(i))),P(g(l(o(i))),g(r(o(i))))),(P(g(i),g(o(i)))))); }; }; }; };
  
void splitTriangles() {    // splits each tirangle into 4
  for (int i=0; i<3*nt; i=i+3) {
    V[3*nt+i]=v(i); V[n(3*nt+i)]=w(p(i)); V[p(3*nt+i)]=w(n(i));
    V[6*nt+i]=v(n(i)); V[n(6*nt+i)]=w(i); V[p(6*nt+i)]=w(p(i));
    V[9*nt+i]=v(p(i)); V[n(9*nt+i)]=w(n(i)); V[p(9*nt+i)]=w(i);
    V[i]=w(i); V[n(i)]=w(n(i)); V[p(i)]=w(p(i));
    };
  nt=4*nt; nc=3*nt;  };
  
void refine() { updateON(); splitEdges(); bulge(); splitTriangles(); updateON();}
  
//  ========================================================== FILL HOLES ===========================================
void fanHoles() {for (int cc=0; cc<nc; cc++) if (visible[t(cc)]&&b(cc)) fanThisHole(cc); normals();  }
void fanThisHole() {fanThisHole(cc);}
void fanThisHole(int cc) {   // fill shole with triangle fan (around average of parallelogram predictors). Must then call computeO to restore O table
 if(!b(cc)) return ; // stop if cc is not facing a border
 G[nv].set(0,0,0);   // tip vertex of fan
 int o=0;              // tip corner of new fan triangle
 int n=0;              // triangle count in fan
 int a=n(cc);          // corner running along the border
 while (n(a)!=cc) {    // walk around the border loop 
   if(b(p(a))) {       // when a is at the left-end of a border edge
      G[nv].add( P(P(g(a),g(n(a))),P(g(a),V(g(p(a)),g(n(a))))) ); // add parallelogram prediction and mid-edge point
      o=3*nt; V[o]=nv; V[n(o)]=v(n(a)); V[p(o)]=v(a); visible[nt]=true; nt++; // add triangle to V table, make it visible
      O[o]=p(a); O[p(a)]=o;        // link opposites for tip corner
      O[n(o)]=-1; O[p(o)]=-1;
      n++;}; // increase triangle-count in fan
    a=s(a);} // next corner along border
 G[nv].mul(1./n); // divide fan tip to make it the average of all predictions
 a=o(cc);       // reset a to walk around the fan again and set up O
 int l=n(a);   // keep track of previous
 int i=0; 
 while(i<n) {a=s(a); if(v(a)==nv) { i++; O[p(a)]=l; O[l]=p(a); l=n(a);}; };  // set O around the fan
 nv++;  nc=3*nt;  // update vertex count and corner count
 };

// =========================================== GEODESIC MEASURES, DISTANCES =============================
void computeDistance(int maxr) { // marks vertices and triangles at a graph distance of maxr
  int tc=0; // triangle counter
  int r=1; // ring counter
  for(int i=0; i<nt; i++) {Mt[i]=0;};  // unmark all triangles
  Mt[t(cc)]=1; tc++;                   // mark t(cc)
  for(int i=0; i<nv; i++) {Mv[i]=0;};  // unmark all vertices
  while ((tc<nt)&&(r<=maxr)) {  // while not finished
     for(int i=0; i<nc; i++) {if ((Mv[v(i)]==0)&&(Mt[t(i)]==r)) {Mv[v(i)]=r;};};  // mark vertices of last marked triangles
     for(int i=0; i<nc; i++) {if ((Mt[t(i)]==0)&&(Mv[v(i)]==r)) {Mt[t(i)]=r+1; tc++;};}; // mark triangles incident on last marked vertices
     r++; // increment ring counter
     };
  rings=r; // sets ring variable for rendring?
  }
  
void computeIsolation() {
  println("Starting isolation computation for "+nt+" triangles");
  for(int i=0; i<nt; i++) {SMt[i]=0;}; 
  for(int c=0; c<nc; c+=3) {println("  triangle "+t(c)+"/"+nt); computeDistance(1000); for(int j=0; j<nt; j++) {SMt[j]+=Mt[j];}; };
  int L=SMt[0], H=SMt[0];  for(int i=0; i<nt; i++) { H=max(H,SMt[i]); L=min(L,SMt[i]);}; if (H==L) {H++;};
  cc=0; for(int i=0; i<nt; i++) {Mt[i]=(SMt[i]-L)*255/(H-L); if(Mt[i]>Mt[t(cc)]) {cc=3*i;};}; rings=255;
  for(int i=0; i<nv; i++) {Mv[i]=0;};  for(int i=0; i<nc; i++) {Mv[v(i)]=max(Mv[v(i)],Mt[t(i)]);};
  println("finished isolation");
  }
  
void clearMt(){
  for(int i=0; i<nt; i++) {Mt[i]=0;}; // reset marking
}

void setMt(int[] tempMt){
  for(int i=0; i<nt;i++){
    Mt[i]=tempMt[i];
  }
}

void addToMt(int[] tempMt){
  for(int i=0; i<nt;i++){
    Mt[i]=max(Mt[i],tempMt[i]);
  }
}

int[] computePath(int startCorner, int endCorner) {                 // graph based shortest path between t(c0 and t(prevc), prevc is the previously picekd corner
  int[] tempMt = new int[maxnt];                 // triangle markers for distance and other things   
  for(int i=0; i<nt; i++) {tempMt[i]=0;}; // reset marking
  tempMt[t(startCorner)]=1; // tempMt[0]=1;            // mark seed triangle

  if(t(startCorner)==t(endCorner)){
    return tempMt;
  }

  for(int i=0; i<nc; i++) {P[i]=false;}; // reset corners as not visited
  int r=1;
  boolean searching=true;
  while (searching) {
     for(int i=0; i<nc; i++) {
       if (searching&&(tempMt[t(i)]==0)&&(!b(i))) { // t(i) is an unvisited triangle and i is not facing a border edge
         if(tempMt[t(o(i))]==r) { // if opposite triangle is ring r
           tempMt[t(i)]=r+1; // mark (invade) t(i) as part of ring r+1
           P[i]=true;    // mark corner i as visited
           if(t(i)==t(endCorner)){searching=false;}; // if we reached the end?
           };
         };
       };
     r++;
     };
  for(int i=0; i<nt; i++) {tempMt[i]=0;};  // graph distance between triangle and t(c)
  rings=1;      // track ring number
  int b=endCorner;
  int k=0;
  while (t(b)!=t(startCorner)) { // back track
    rings++;  
    if (P[b]) {b=o(b); } else {if (P[p(b)]) {b=r(b); } else {b=l(b);};}; tempMt[t(b)]=rings; 
  };
  return tempMt;
}

 void  showDistance() {
   noStroke(); 
   for(int t=0; t<nt; t++) 
     if(Mt[t]!=0) {
       fill(ramp(Mt[t],rings));
       showShrunkOffsetT(t,1,1);
     }; 
   noFill(); 
 } 


//  ==========================================================  GARBAGE COLLECTION ===========================================
void clean() {
   excludeInvisibleTriangles();  println("excluded");
   compactVO(); println("compactedVO");
   compactV(); println("compactedV");
   normals(); println("normals");
   computeO();
   resetMarkers();
   identifyBorderVertices();
   }  // removes deleted triangles and unused vertices
   
void excludeInvisibleTriangles () {for (int b=0; b<nc; b++) {if (!visible[t(o(b))]) {O[b]=b;};};}
void compactVO() {  
  int[] U = new int [nc];
  int lc=-1; for (int c=0; c<nc; c++) {if (visible[t(c)]) {U[c]=++lc; }; };
  for (int c=0; c<nc; c++) {if (!b(c)) {O[c]=U[o(c)];} else {O[c]=c;}; };
  int lt=0;
  for (int t=0; t<nt; t++) {
    if (visible[t]) {
      V[3*lt]=V[3*t]; V[3*lt+1]=V[3*t+1]; V[3*lt+2]=V[3*t+2]; 
      O[3*lt]=O[3*t]; O[3*lt+1]=O[3*t+1]; O[3*lt+2]=O[3*t+2]; 
      visible[lt]=true; 
      lt++;
      };
    };
  nt=lt; nc=3*nt;    
  println("      ...  NOW: nv="+nv +", nt="+nt +", nc="+nc );
  }

void compactV() {  
  println("COMPACT VERTICES: nv="+nv +", nt="+nt +", nc="+nc );
  int[] U = new int [nv];
  boolean[] deleted = new boolean [nv];
  for (int v=0; v<nv; v++) {deleted[v]=true;};
  for (int c=0; c<nc; c++) {deleted[v(c)]=false;};
  int lv=-1; for (int v=0; v<nv; v++) {if (!deleted[v]) {U[v]=++lv; }; };
  for (int c=0; c<nc; c++) {V[c]=U[v(c)]; };
  lv=0;
  for (int v=0; v<nv; v++) {
    if (!deleted[v]) {G[lv].set(G[v]);  deleted[lv]=false; 
      lv++;
      };
    };
 nv=lv;
 println("      ...  NOW: nv="+nv +", nt="+nt +", nc="+nc );
  }

// ============================================================= ARCHIVAL ============================================================
boolean flipOrientation=false;            // if set, save will flip all triangles

void saveMesh() {
  String [] inppts = new String [nv+1+nt+1];
  int s=0;
  inppts[s++]=str(nv);
  for (int i=0; i<nv; i++) {inppts[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z);};
  inppts[s++]=str(nt);
  if (flipOrientation) {for (int i=0; i<nt; i++) {inppts[s++]=str(V[3*i])+","+str(V[3*i+2])+","+str(V[3*i+1]);};}
    else {for (int i=0; i<nt; i++) {inppts[s++]=str(V[3*i])+","+str(V[3*i+1])+","+str(V[3*i+2]);};};
  saveStrings("mesh.vts",inppts);  println("saved on file");
  };

void loadMesh() {
  println("loading fn["+fni+"]: "+fn[fni]); 
  String [] ss = loadStrings(fn[fni]);
  String subpts;
  int s=0;   int comma1, comma2;   float x, y, z;   int a, b, c;
  nv = int(ss[s++]);
    print("nv="+nv);
    for(int k=0; k<nv; k++) {int i=k+s; 
      comma1=ss[i].indexOf(',');   
      x=float(ss[i].substring(0, comma1));
      String rest = ss[i].substring(comma1+1, ss[i].length());
      comma2=rest.indexOf(',');    y=float(rest.substring(0, comma2)); z=float(rest.substring(comma2+1, rest.length()));
      G[k].set(x,y,z);
    };
  s=nv+1;
  nt = int(ss[s]); nc=3*nt;
  println(", nt="+nt);
  s++;
  for(int k=0; k<nt; k++) {int i=k+s;
      comma1=ss[i].indexOf(',');   a=int(ss[i].substring(0, comma1));  
      String rest = ss[i].substring(comma1+1, ss[i].length()); comma2=rest.indexOf(',');  
      b=int(rest.substring(0, comma2)); c=int(rest.substring(comma2+1, rest.length()));
      V[3*k]=a;  V[3*k+1]=b;  V[3*k+2]=c;
    }
  }; 
  

void loadMeshOBJ() {
  println("loading heart.obj"); 
  String [] ss = loadStrings("heart.obj");
  String subpts;
  String S;
  int comma1, comma2;   float x, y, z;   int a, b, c;
  int s=2;   
  println(ss[s]);
  int nn=ss[s].indexOf(':')+2; println("nn="+nn);
  nv = int(ss[s++].substring(nn));  println("nv="+nv);
  int k0=s;
    for(int k=0; k<nv; k++) {int i=k+k0; 
      S=ss[i].substring(2); if(k==0 || k==nv-1) println(S);
      comma1=S.indexOf(' ');   
      x=float(S.substring(0, comma1));
      String rest = S.substring(comma1+1);
      comma2=rest.indexOf(' ');    y=float(rest.substring(0, comma2)); z=float(rest.substring(comma2+1));
      G[k].set(x,y,z); if(k<3 || k>nv-4) {print("k="+k+" : "); }
      s++;
    };
  s=s+2; 
  println("Triangles");
  println(ss[s]);
  nn=ss[s].indexOf(':')+2;
  nt = int(ss[s].substring(nn)); nc=3*nt;
  println(", nt="+nt);
  s++;
  k0=s;
  for(int k=0; k<nt; k++) {int i=k+k0;
      S=ss[i].substring(2);                        if(k==0 || k==nt-1) println(S);
      comma1=S.indexOf(' ');   a=int(S.substring(0, comma1));  
      String rest = S.substring(comma1+1); comma2=rest.indexOf(' ');  
      b=int(rest.substring(0, comma2)); c=int(rest.substring(comma2+1));
      V[3*k]=a-1;  V[3*k+1]=b-1;  V[3*k+2]=c-1;
    }
  for (int i=0; i<nv; i++) G[i].mul(4);  
  }; 

 
  } // ==== END OF MESH CLASS
  
vec labelD=new vec(-4,+4, 12);           // offset vector for drawing labels  

