
/* output */
void Output_gnuplot( FILE *gp, double *X, double point_size, double arrowlength, double ptime, double v0, double u0)
{
  double *x = X, *y = X+N, *theta = X+4*N;
  fprintf(gp, "set title 'N = %d, time = %.2f v0 = %.2f u0 = %.2f'\n",N,ptime,v0,u0);
  //fprintf(gp, "set style fill\n");
  //fprintf(gp, "plot '-' pt 7 ps %f lt 7, '-' with lines lt -1\n",point_size);
  fprintf(gp, "plot '-' with ellipse lt 7, '-' with lines lt -1\n");
  for(int i=0; i<N; i++){
    fprintf(gp,"%f\t%f\t%f\n", x[i], y[i],point_size);    // データの書き込み
  }
  fprintf(gp,"e\n");
  for(int i=0; i<N; i++){
    fprintf(gp,"%f\t%f\n%f\t%f\n\n", x[i], y[i],x[i]+arrowlength*cos(theta[i]),y[i]+arrowlength*sin(theta[i]));    // データの書き込み
  }
  fprintf(gp,"e\n");
  fflush(gp);
  usleep(600000);
}

