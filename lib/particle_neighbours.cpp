
class Neighbourlist{
  public:
    int  i,j;
    char bx, by;
    Neighbourlist(int ci, int cj, char cbx, char cby)
       {  i= ci; j = cj; bx = cbx; by = cby;  }
};
std::vector<Neighbourlist> neighbours;


void ReConfiguration(double *x,double *y)
{
  /**  **/
  for(int i=0; i<N; ++i)x[i] = fmod(x[i]+LX,LX);
  for(int i=0; i<N; ++i)y[i] = fmod(y[i]+LY,LY);
  /**
  for(int i=0; i<N; ++i){
    if( x[i]<0 )x[i]+=LX;
    else if(x[i]>=LX) x[i]-=LX;
  }
  for(int i=0; i<N; ++i){
    if( y[i]<0 )y[i]+=LY;
    else if(y[i]>=LY) y[i]-=LY;
  }
  **/
}



void makeneighbours(double *x,double *y)
{
  ReConfiguration(x,y);

  // make cell_list
  std::vector<int> cell_list[cell_size_x*cell_size_y];
  for( int i=0; i<N; ++i ){
    cell_list[(int)(x[i]/h)*cell_size_y+(int)(y[i]/h)].push_back(i);
  }

  // make Verlet list
  std::vector<Neighbourlist> neighboursTemp[cell_size_x];
  for( int I=0; I<cell_size_x; ++I ){
    neighboursTemp[I].reserve(5*N);

    for( int J=0; J<cell_size_y; ++J ){
      int oI = I*cell_size_y + J;
      // self
      for( auto itr=cell_list[oI].begin(); itr!=cell_list[oI].end(); ++itr ) {
        for( auto ptr=itr+1; ptr!=cell_list[oI].end(); ++ptr ){
          double dx = x[*ptr]-x[*itr];
          double dy = y[*ptr]-y[*itr];
          if( dx*dx+dy*dy < verlet_radius_sq )neighboursTemp[I].emplace_back(*itr,*ptr,0,0);
        }
      }

      int tI2 = (I+1-(I==cell_size_x-1)*cell_size_x)*cell_size_y + J;
      for(const auto& e0 : cell_list[oI] ){
        for(const auto& e1 : cell_list[tI2] ){
          double dx = x[e1]-x[e0] + LX*(I==cell_size_x-1);
          double dy = y[e1]-y[e0];
          if( dx*dx+dy*dy < verlet_radius_sq)neighboursTemp[I].emplace_back(e0,e1,+(I==cell_size_x-1),0);
        }
      }

      int tI3 = (I+1-(I==cell_size_x-1)*cell_size_x)*cell_size_y + J+1-(J==cell_size_y-1)*cell_size_y;
      for(const auto& e0 : cell_list[oI] ){
        for(const auto& e1 : cell_list[tI3] ){
          double dx = x[e1]-x[e0] + LX*(I==cell_size_x-1);
          double dy = y[e1]-y[e0] + LY*(J==cell_size_y-1);
          if( dx*dx+dy*dy < verlet_radius_sq)neighboursTemp[I].emplace_back(e0,e1,+(I==cell_size_x-1),+(J==cell_size_y-1));
        }
      }

      int tI4 = I*cell_size_y + J+1-(J==cell_size_y-1)*cell_size_y;
      for(const auto& e0 : cell_list[oI] ){
        for(const auto& e1 : cell_list[tI4] ){
          double dx = x[e1]-x[e0];
          double dy = y[e1]-y[e0] + LY*(J==cell_size_y-1);
          if( dx*dx+dy*dy < verlet_radius_sq)neighboursTemp[I].emplace_back(e0,e1,0,+(J==cell_size_y-1));
        }
      }

    }
  }

  neighbours.clear();
  int totN=0;
  for(int I=0; I<cell_size_x; I++)totN += (int)neighboursTemp[I].size();
  neighbours.reserve( totN );
  for(int I=0; I<cell_size_x; I++){
      std::copy(neighboursTemp[I].begin(),neighboursTemp[I].end(),std::back_inserter(neighbours));
  }
}



extern double MARGIN;
inline void check_pairlist(double vx[],double vy[],double x[],double y[] )
{
  static double margin_length = MARGIN;
  double maxV = 0.0;
  for(int i=0; i<N; i++ ){
    if(maxV < vx[i]*vx[i]+vy[i]*vy[i] )maxV=vx[i]*vx[i]+vy[i]*vy[i];
  }
  maxV = sqrt(maxV);
  margin_length -= maxV*2.0*Dt;
  if(margin_length < 0.0){
    margin_length = MARGIN;
    makeneighbours(x,y);
  }
}


