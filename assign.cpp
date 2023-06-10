

/**
 * \mainpage Decomposition of Polygon to Convex Polygons
 *
 *   \section Soutik Soutik          - 2019B5A71392H
     \section Sahil Sahil           - 2019B5A71379H
     \section Tarun Tarun           - 2019B3A70733H
     \section Shashank Shashank   - 2019B4A70956H
            
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include<bits/stdc++.h>

using namespace std;
 
struct Vertex;
struct Edge;
//struct Face;
FILE* fpi=fopen("output.txt","w");
struct Vertex {                 ///Struct to store a vertex
    double x, y;                ///x and y co-ordinates of a vertex
    Edge* inc_edge1=nullptr;    ///the incident edge that runs clock-wise
    Edge* inc_edge2=nullptr;    ///the incident edge that runs anti-clockwise
    
    Vertex (){
        x=0.0;
        y=0.0;
        inc_edge1=nullptr;
        inc_edge2=nullptr;
    }
    
    Vertex(pair<double,double> a){      ///Makes a vertex with the given co-ordinates
        x=a.first;
        y=a.second;
        inc_edge1 = nullptr;   
        inc_edge2 = nullptr;    
    }
};


struct LUP{                     ///Struct used to merge polygons with the list of diagonals and the corresponding sub-polygons it belongs to
    int i;
    Vertex* v1;
    Vertex* v2;
};


double angle(Vertex* , Vertex* , Vertex* );
struct Edge {                   ///Struct to store an Edge
    Vertex* origin=nullptr;     ///Origin Vertex of the edge
    Vertex* destination=nullptr;///Destination vertex of the edge
    Edge* twin=nullptr;         ///Twin edge corresponding to this edge
    // Edge* next = nullptr;
    // Edge* prev = nullptr;
    // Face* face = nullptr;
    Edge(){
       origin=nullptr;
       destination=nullptr;
       twin=nullptr;
    }
    
    Edge(Vertex* a,Vertex* b){  ///Makes an Edge between given two vertices
        origin= a;
        destination =b;
        Edge* twin=nullptr;
    }
};

// struct Face {
//     Edge* edge = nullptr;
// };

vector<pair<double,double>> make_list(vector<Vertex*> );

class DCEL {                    ///Class which will contain the polygon
public:
    vector<Vertex*> vertices;   ///Vector containing all of the vertices in the polygon
//    vector<Edge*> edges;
//    vector<Face*> faces;
    
    DCEL (vector<pair<double,double>> vert){    ///Constructor that makes a polygon out of the list of co-ordinates given
        for(pair<double,double> a:vert){
            Vertex* temp=new Vertex(a);         ///Making a new vertex
            vertices.push_back(temp);
        }
        int length=vertices.size();
        for(int i=0;i<length;i++){              ///Making Edges between all of the vertices
            Edge* temp1=new Edge(vertices[i],vertices[(i+1)%length]);
            Edge* temp2=new Edge(vertices[(i+1)%length],vertices[i]);
            temp1->twin=temp2;
            temp2->twin=temp1;
            vertices[(i+1)%length]->inc_edge1=temp1;
            vertices[i]->inc_edge2=temp2;
            //edges.push_back(temp1);
            //edges.push_back(temp2);
        } 
    }
    
    bool notches (vector<pair<double,double>> vect){    ///Finds out whether the calling polygon has a notch inside it
        double min_x,max_x=vect[0].first;               ///Creating the rectangle 
        double min_y,max_y=vect[0].second;
        for(pair<double,double> vert:vect){
            if(max_x<vert.first)max_x=vert.first;
            if(max_y<vert.second)max_y=vert.second;
            if(min_x>vert.first)min_x=vert.first;
            if(min_y>vert.second)min_y=vert.second;
        }
        
        for(int i=0;i<(int)vect.size();i++){
            if(vect[i].first<min_x||vect[i].first>max_x||vect[i].second<min_y||vect[i].second>max_y)continue;   ///Checking whether the vertex lies inside the rectangle or not
            int intersections=0;                                                                                ///The vertex is found out to be lying inside the rectangle, so drawing a horizontal line to the right starting from it, and counting the number of intersections with the edges of the polygon
            int n=vertices.size();
            for (int j = 0; j < n; j++) {
                Vertex* p1 = vertices[j];
                Vertex* p2 = vertices[(j + 1) % n];

                if ((vect[i].second >= min(p1->y, p2->y)) && (vect[i].second <= max(p1->y, p2->y))) {
                    if (vect[i].first <= max(p1->x, p2->x)) {
                        if (p1->y != p2->y) {
                            double xIntersection = (vect[i].second - p1->y) * (p2->x - p1->x) / (p2->y - p1->y) + p1->x;
                            if (p1->x == p2->x || vect[i].first <= xIntersection) {
                                intersections++;
                            }
                        }
                    }
                }
            }
            if(intersections%2){remove_notch(vect[i]);return false;}            ///If there are an odd number of intersections, then the vertex is a notch, otherwise not
        }
        return true;
    }
    
    void remove_notch(pair<double,double> p1){          ///Removes the notch from the calling polygon
        
        vector<Vertex*> rem;
        int n=vertices.size();
        
        Vertex* p2=vertices[0];
        Vertex* p3=vertices[n-1];
        
        rem.push_back(p3);
        double value=(p3->y-p1.second)-(p3->x-p1.first)*(p2->y-p1.second)/(p2->x-p1.first);
        
        for(Vertex* p3:vertices){
            double new_val=(p3->y-p1.second)-(p3->x-p1.first)*(p2->y-p1.second)/(p2->x-p1.first);   ///Checking whether any vertex lies on the same side as the last added vertex of the line between the first added vertex and the notch 
            if(value*new_val>0)rem.push_back(p3);                                                   ///If a vertex is on the same side as the last added vertex, pushing it into a vector
        }
        
        for(Vertex* p3:rem)remove_vertex(make_pair(p3->x,p3->y));                                   ///Removing all the vertices existing in the aforementioned vector from the polygon
        
    }
    
    void remove_vertex(pair<double,double> p){  ///removes the vertex that has the same co-ordinates as the pair in the parameter
        int n=vertices.size();
        Vertex* v;
        Vertex* v1;
        Vertex* v2;
        for(int i=0;i<n;i++){
            v=vertices[i];
            if(v->x==p.first&&v->y==p.second){
                v1 = (v->inc_edge1)->origin;
                v2 = (v->inc_edge2)->origin;
                vertices.erase(vertices.begin()+i); ///Removing the vertex from the list
                break;
            }
        }
        Edge* Edge1=new Edge(v1,v2);
        Edge* Edge2=new Edge(v2,v1);
        // delete (v->inc_edge1);
        // delete (v->inc_edge2);
        // delete (v1->inc_edge2);
        // delete (v2->inc_edge1);
        // delete v;
        v2->inc_edge1=Edge1;                    ///Making new edges, and over-writing the edges of the vertices on each side with the newly created ones
        v1->inc_edge2=Edge2;
        Edge1->twin=Edge2;
        Edge2->twin=Edge1;
    }
    
    void add_vertex(pair<double,double> p){     ///adds a vertex at the last position for a polygon
        Vertex* new_vertex=new Vertex(p);       ///making a new vertex
        vertices.push_back(new_vertex);
        Vertex* v1=vertices[(int)vertices.size()-2];
        Vertex* v2=vertices[0];
        Edge* e1=new Edge(v1,new_vertex);       ///Making new edges and setting pointers accordingly
        Edge* e2=new Edge(new_vertex,v1);
        Edge* e3=new Edge(v2,new_vertex);
        Edge* e4=new Edge(new_vertex,v2);
        v1->inc_edge2=e2;
        new_vertex->inc_edge1=e1;
        v2->inc_edge1=e4;
        new_vertex->inc_edge2=e3;
        //delete (v1->inc_edge2);
        //delete (v2->inc_edge1);
    } 

    void print_polygon(){               ///prints the polygon(testing purposes)
        printf("Vertices vector\n");
        for(Vertex* v:vertices)printf("%lf %lf  ",v->x,v->y);
        printf("\nClockwise polygon\n");
        Vertex* v=vertices[0];
        Vertex* v1=vertices[0];
        do{
            printf("%lf %lf -> ",v->x,v->y);
            v=v->inc_edge2->origin;
        }while(v1!=v);
        printf("\nCounterclockwise\n");
        v=vertices[0];
        v1=vertices[0];
        do{
            printf("%lf %lf -> ",v->x,v->y);
            v=v->inc_edge1->origin;
        }while(v1!=v);
    }
    void print_polygon_file(){               ///prints the polygon(in input.txt)
        for(Vertex* v:vertices)fprintf(fpi,"%lf %lf  ",v->x,v->y);
        fprintf(fpi,"\n");
    }

    pair<Vertex*,Vertex*> remove_polygon(DCEL sub_polygon){         ///Removes the polygon in the parameter from the calling polygon(calling polygon)
        
        Vertex *v1;
        Vertex *v2;
        int n = sub_polygon.vertices.size();
        int n1=vertices.size();
        int j=0;
        Vertex* vert1=vertices[0];
        while(vert1->x!=sub_polygon.vertices[0]->x||vert1->y!=sub_polygon.vertices[0]->y)vert1=vertices[++j];       ///Finding the vertices between which a new diagonal was formed
        
        int i;
        for(i=j;i<n1;i++){
            v1=vertices[i%n1];
            v2=sub_polygon.vertices[(i-j)%n];
            if(v1->x!=v2->x||v1->y!=v2->y)break;
        }

        int x1=v2->x;
        int x2=v2->inc_edge1->origin->x;
        int y1=v2->y;
        int y2=v2->inc_edge1->origin->y;

        int flag=0;

        for(i=0;i<n1;i++){
            if(flag==2)break;
            if((vertices[i]->x==x1&&vertices[i]->y==y1)||(vertices[i]->x==x2&&vertices[i]->y==y2)){
                if(flag==0)v1=vertices[i],flag=1;
                else v2=vertices[i],flag=2;
            }
        }

        vector<pair<double,double>> sub=make_list(sub_polygon.vertices);           

        for(i=0;i<(int)sub.size();i++){             ///Removing the vertices containing the diagonal from the sub-polygon
            if(((sub[i].first==v1->x)&&(sub[i].second==v1->y))||((sub[i].first==v2->x)&&(sub[i].second==v2->y))){
               sub.erase(sub.begin()+i);
            }
        }

        for(pair<double,double> rem:sub){
            for(i=0;i<n1;i++)if(vertices[i]->x==rem.first&&vertices[i]->y==rem.second)vertices.erase(vertices.begin()+i);   ///Removing all the remianing vertices in the sub-polygon from the original polygon
        }

        Edge* temp1=new Edge(v1,v2);            ///Making and adding new edges to our vertices
        Edge* temp2=new Edge(v2,v1);
        temp1->twin=temp2;
        temp2->twin=temp1;
        v2->inc_edge1=temp1;
        v1->inc_edge2=temp2;

        return make_pair(v2,v1);
    }

    bool convexpolygon(){                       ///Checks whether calling polygon is convex or not
        Vertex* a=vertices[0];
        Vertex* b=vertices[1];
        Vertex* c=vertices[2];
        int i=2;
        int j=2;
        bool conv=true;
        while(j<=vertices.size()+1){            ///Checks all inside angles for convexity
            if(angle(a,b,c)<M_PI){
                a=b;b=c;c=vertices[(++i)%vertices.size()];
            }
            else {conv=false;break;}
            j++;
        }
        return conv;
    }

};

double angle(Vertex* a, Vertex* b, Vertex* c) {         ///Finds and returns the counter clockwise angle (in radian) between the lines ab and bc
    double dx1 = a->x - b->x;
    double dy1 = a->y - b->y;
    double dx2 = c->x - b->x;
    double dy2 = c->y - b->y;
    double dot_product = dx1 * dx2 + dy1 * dy2;
    double cross_product = dx1 * dy2 - dx2 * dy1;
    double angle = atan2(cross_product, dot_product);
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return angle;
}

DCEL make_polygon(vector<Vertex*> list){    ///Takes a Vector containing Vertices and converts it into a polygon.
                                            ///Used to convert the list generated by MP1 to a polygon
    vector<pair<double,double>> poly;

    for(int i=0;i<list.size();i++){
        poly.push_back(pair<double,double>(list[i]->x,list[i]->y));
    }

    DCEL polygon(poly);

    return polygon;
}

vector<pair<double,double>> make_list(vector<Vertex*> list){    ///Converts the vector containing Vertices to a vector containing Coordinates of the polygon
                                                                ///Used to remove notches in MP1.
    vector<pair<double,double>>  coord;

    for(int i=0;i<list.size();i++){
        coord.push_back(make_pair(list[i]->x,list[i]->y));
    }

    return coord;
}

vector<Vertex*> MP1(DCEL polygon){          ///Start of the MP1 procedure as defined in the paper
    
    int n = polygon.vertices.size();

    if(n<4){
        return polygon.vertices;            ///If number of vertices is less than 4 it must be a triangle to be a polygon which can't be concave
    }

    vector<Vertex*> list;
    vector<pair<double,double>> notch;

    list.push_back(polygon.vertices[0]);    
    list.push_back(polygon.vertices[1]);

    Vertex* v1 = list[0];
    Vertex* v2 = list[1];
    Vertex* vim1 = v1;
    Vertex* vi = v2;
    Vertex* vip1 = polygon.vertices[2];
    int i=3;

    while(list.size()<polygon.vertices.size()){
        if(angle(vim1,vi,vip1)<M_PI && angle(vi,vip1,v1)<M_PI && angle(vip1,v1,v2)<M_PI){        ///Check if the angle to be added to the list of convex polygon is Convex by checking with all three angles that can be formed with the added vertex.
            list.push_back(vip1);                                                                   ///If not then we stop the algorithm and return the Polygon formed
        }
        else{
            break;
        }

        vim1=vi;
        vi=vip1;
        vip1 = polygon.vertices[((i++)%n)];
    }

    if(n == list.size()){                           ///If the list generated is the same size as the polygon passed then there are no concave angles and we return the list without further checking.
        return list;
    }
    else if(list.size() == 2){                      ///If the list generated is only of size 2 then it's not a polygon and we return the list which just resends the polygon with new edges.
        cout<<"Edge\n";
        return list;
    }

    else{

        DCEL TEMP = make_polygon(list);

        vector<pair<double,double>> poly_list,sub_poly_list;
        poly_list = make_list(polygon.vertices);
        sub_poly_list = make_list(list);    
                                                                            ///Removing vertices from the original polygon to get the vertices not present in the generated list
        poly_list.erase( remove_if( begin(poly_list),end(poly_list),[&](auto x){return find(begin(sub_poly_list),end(sub_poly_list),x)!=end(sub_poly_list);}), end(poly_list) );
                                                                                                                
        while(!TEMP.notches(poly_list)){        ///Checking to find notches in the newly generated sub-polygon using the vertices not present in the polygon to check
            cout<<"Removing notches";
        }
        
        return TEMP.vertices;                   ///Returning the polygon after removing the notches
    }
}

bool compare(pair<double,double> l1,pair<double,double> l2){
    if(l1.first==l2.first && l1.second==l2.second){
        return 1;
    }    
    else{
        return 0;
    }
}

DCEL Merger(DCEL p1, DCEL p2){///Merger function takes in graphs as input
        vector<pair<double,double>> final;
        vector<pair<double,double>> final2;
        vector<pair<double,double>> first=make_list(p1.vertices);
        vector<pair<double,double>> second=make_list(p2.vertices);
        final.reserve( first.size() + second.size() ); 
        final.insert( final.end(), first.begin(), first.end() );
        final.insert( final.end(), second.begin(), second.end() );///merging the list of vertices of both subgraphs p1 and p2
            for(int it=0;it<final.size();it++){
               for(int it2=0;it2!=final.size();it2++){
                   if(it!=it2){
                   if(compare(final[it],final[it2])){
                       final.erase(final.begin()+it2);
                       it2--;
                   }
               } }
            }/// removing the duplicate elements (which would be the diangonal)
        DCEL p(final);///creating new graph p by combining p1 and p2
        return p;///output of merger is the new merged graph

}

vector<DCEL> merging(vector<pair<Vertex*,Vertex*>> LLE,vector<LUP> LP,DCEL p,vector<DCEL> subpolygons){/// merging function takes in the list of diagonals(LLE),a list of polygons containing the corresponding diagonal(LP), the orignal polygon(P) and list of subpolygons made by process mp1(subpolygons)
    int m=LLE.size();//no. of elements in lle
    int NP=m+1;//no. of polygons;
    vector<pair<int,int>> LDP;//all true initially
    for(int j=0;j<m;j++){///for all elements in LLE structure,i.e, for all diagonals
        Vertex* vs;
        Vertex* vt;
        vs=LLE[j].first;
        vt=LLE[j].second;
        int k1=0,k2=0;///k1 and k2 are the indices of subpolygons in subpolygons list,which correspond to the polygonson either side of diagonal vs,vt 
        
        
        for(auto it=LP.begin();it!=LP.end();++it){
                if(it->v1==vs && it->v2==vt){
                    k1=it->i;}
                else if(it->v1==vt && it->v2==vs){
                    k2=it->i;
                    }
                }
        

        DCEL F=Merger(subpolygons[k1],subpolygons[k2]);///F is the merged polygon after removing the diagonal between K1 and K2,this is done by the Merger function
        if(F.convexpolygon()){  ///checking if the new merged polygon(F) is convex or not      
                    NP=NP+1;
                    pair<int,int> p=make_pair(k1,k2);
                    LDP.push_back(p);
                    //int s1=k1;
                    //subpolygons[s1]=F;
                    //auto it=subpolygons.begin()+k2;
                    //subpolygons.erase(it);
                         
        }
    }
    for(auto it=LDP.begin();it!=LDP.end();++it){
        int s1=it->first;
        subpolygons[s1]=Merger(subpolygons[it->first],subpolygons[it->second]);
        
        ///adding  merged polygon at position of p1 to subpolygon list and remove p1,p2 
    }
    for(auto it=LDP.begin();it!=LDP.end();++it){
        auto it2=subpolygons.begin()+(it->second);
        subpolygons.erase(it2);
    }

    return subpolygons;///merging function returns the list of new subpolygons after merging

}

int main(){

    int n;
    printf("Number of vertices\n");
    cin>>n;
    vector<pair<double,double>> vertex(n); 
    cout<<"Enter vertices in 'x y' format ,where x and y are coordinates of the vertice,enter them in clockwise order"<<endl;
    for(int i=0;i<n;i++){
        double a,b;
        cin>>a;
        cin>>b;
        vertex[i]=make_pair(a,b);
    }
    
    DCEL polygon(vertex);
    polygon.print_polygon();


    DCEL polygon_1(make_list(polygon.vertices));
    DCEL sub_poly(make_list(polygon.vertices));
    vector<Vertex*> temp;
    vector<DCEL> sub_polygons;
    pair<Vertex*,Vertex*> diagonals,prev_diagonal;
    vector<LUP> merge_help; 
    vector<pair<Vertex*,Vertex*>> diagonal_list,di_list;
    LUP lup_last;
    
    do{

        temp=MP1(polygon_1);                ///Calling the procedure MP1

        if(temp.size()==polygon_1.vertices.size()){
            sub_polygons.push_back(polygon_1);          ///If the subpolygon is the same size as the input polygon it must be the final decomposition.
            break;
        }

        else if(temp.size()>2){

            sub_poly=make_polygon(temp);
            sub_polygons.push_back(sub_poly);           ///Removing the Sub-Polygons from the Polygon to recursively generate decompositions 
            diagonals=polygon_1.remove_polygon(sub_poly);

            diagonal_list.push_back(make_pair(diagonals.first,diagonals.second));       ///Storing the diagonals generated to use in the Merging phase
            di_list.push_back(make_pair(diagonals.first,diagonals.second));
            di_list.push_back(make_pair(diagonals.second,diagonals.first));

        }

        else{
            ///Pop the first vertex of polygon_1 and add it to the end to reorder the list of vertices given to MP1.
            polygon_1.vertices.push_back(polygon_1.vertices[0]);
            polygon_1.vertices.erase(polygon_1.vertices.begin());
        }    

    }while(polygon_1.vertices.size()>0);

    // cout<<"diagonal_list\n";
    // for(int i=0;i<diagonal_list.size();i++){
    //     cout<<diagonal_list[i].first->x<<" "<<diagonal_list[i].first->y<<" "<<diagonal_list[i].second->x<<" "<<diagonal_list[i].second->y;
    //     cout<<"\n";
    // }
    // cout<<"the lup list\n";
    
    for(int i=0;i<di_list.size();i++){
        for(int j=0;j<sub_polygons.size();j++){                 ///Geneating the list LUP as described in the paper
            for(int k=0;k<sub_polygons[j].vertices.size();k++){
                if((sub_polygons[j].vertices[k]->x==di_list[i].first->x) && (sub_polygons[j].vertices[k]->y==di_list[i].first->y) && 
                (sub_polygons[j].vertices[(k+1)%(sub_polygons[j].vertices.size())]->x==di_list[i].second->x) && 
                (sub_polygons[j].vertices[(k+1)%(sub_polygons[j].vertices.size())]->y==di_list[i].second->y)){
                    LUP lup = {j,di_list[i].first,di_list[i].second};
                    merge_help.push_back(lup);
                    
                }
            }
        }
    }
    
    // for(int i=0;i<merge_help.size();i++){
    //     cout<<merge_help[i].i<<" "<<merge_help[i].v1->x<<" "<<merge_help[i].v1->y<<" "<<merge_help[i].v2->x<<" "<<merge_help[i].v2->y<<"\n";
    // }

    cout<<"\nOUTPUT BEFORE MERGING\n";

    for(int i=0;i<sub_polygons.size();i++){
        sub_polygons[i].print_polygon();
        cout<<"\n\n";
    }

    ///Calling the Merging process
    sub_polygons=merging(diagonal_list,merge_help,polygon,sub_polygons);

    /// Printing Final Polygon
    printf("\n OUTPUT AFTER MERGING PROCESS\n");
    for(int i=0;i<sub_polygons.size();i++){
        sub_polygons[i].print_polygon();
        cout<<"\n\n";
    }


    polygon.print_polygon_file();
    for(int i=0;i<sub_polygons.size();i++){
        sub_polygons[i].print_polygon_file();
    }


    return 0;
}