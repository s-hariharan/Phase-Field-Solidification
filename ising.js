var canvas = document.getElementById('theCanvas');
var context = canvas.getContext('2d');
var image = context.createImageData(canvas.width, canvas.height);




var i,j,k,t,b,l,r;
var nx=100;
var ny=100;
var tstep=100;
            //initialization of non dimensionalized variables.
var delx=0.03;
var dely=0.03;
var alpha=0.9;
var epsilonbar=0.01;
var tau=0.0003;
var delt=0.00002;
var lamda=10;
var pi=3.14;
var m;
var g_phi;
var kappa;
var lheat=2;
var size = nx;
var noise_str = 0.01;
var squareWidth = canvas.width/size;
var delta = 0.005;
var mode = 4;
var counter = 0;


var phi = new Array(size);
var temp = new Array(size);
var new_phi = new Array(size);
var new_temp = new Array(size);
var grad_x = new Array(size);
var grad_y = new Array(size);
var angl = new Array(size);


var startTime = (new Date()).getTime();

var running= false;
for(i=0;i<nx;i++)
{ phi[i] = new Array(size);
  temp[i]= new Array(size);
  new_phi[i] = new Array(size);
  new_temp[i]= new Array(size);
  grad_x[i] = new Array(size);
  grad_y[i] = new Array(size);
  angl[i] = new Array(size);
  for(j=0;j<ny;j++)
  {  temp[i][j]=0;
     angl[i][j]=0;

      if((i-nx/2-1)*(i-nx/2-1)+(j-100)*(j-100)<10)
      {
     phi[i][j]=1;
     temp[i][j]=1;
     }
     else
     {phi[i][j]=0;}
     colorSquare(i,j);
  }
}

context.putImageData(image,0,0);
simulate();

function simulate(){
  if(running){
    evolver();
    for(i=0;i<nx;i++){
      for(j=0;j<ny;j++){
        colorSquare(i,j);
      }
    }
    context.putImageData(image,0,0);

  }
  window.setTimeout(simulate,200);
}


 function colorSquare(i,j) {

  for(py=j*squareWidth; py<(j+1)*squareWidth; py++) {
    for( px= i*squareWidth; px<(i+1)*squareWidth; px++){
      var index = (px +py*image.width)*4;
      image.data[index+0] = Math.floor(phi[i][j]*240);
      image.data[index+1] = 0;
      image.data[index+2] = 125;
      image.data[index+3] = 255;
    }
  }
}


function evolver () {
   for( k =0; k<100;k++){
    for(i=0;i<nx;i++)
    {
       for(j=0;j<ny;j++)
        {
            var angle=lamda*(1-temp[i][j]);
            m=(alpha/pi)*Math.atan(angle);
            g_phi=phi[i][j]*(1-phi[i][j])*(phi[i][j]-0.5+m) ;
            b=j-1;
            t=j+1;
            r=i+1;
            l=i-1;

            if(i==0)
                l=i;
            if(i==nx-1)
                r=i;
            if(j==0)
                b=j;
            if(j==ny-1)
                t=j;

            grad_x[i][j]=(phi[r][j]-phi[l][j])/(2*delx);
            grad_y[i][j]=(phi[i][t]-phi[i][b])/(2*dely);


               if(grad_x[i][j]==0)
               {
                  if(grad_y[i][j]<0)
                     angl[i][j]=-0.5*pi;
                  else if(grad_y[i][j]>0)
                     angl[i][j]=0.5*pi;
                  else
                     angl[i][j]=0;
               }

               if(grad_x[i][j]!=0)
               {
                   angl[i][j]=Math.atan(grad_y[i][j]/grad_x[i][j]);
               }


               var epsilon_rj,epsilon_lj,epsilon_it,epsilon_ib;
               var epsilon_der_rj,epsilon_der_lj,epsilon_der_it,epsilon_der_ib;

            if(phi[i][j]>0|| phi[i][j]<1){
            epsilon_rj=(1+delta*Math.cos(mode*(angl[r][j])))*epsilonbar;
            epsilon_lj=(1+delta*Math.cos(mode*(angl[l][j])))*epsilonbar;
            epsilon_it=(1+delta*Math.cos(mode*(angl[i][t])))*epsilonbar;
            epsilon_ib=(1+delta*Math.cos(mode*(angl[i][b])))*epsilonbar;


            epsilon_der_rj=-1*delta*mode*Math.sin(mode*(angl[r][j]))*epsilonbar;
            epsilon_der_lj=-1*delta*mode*Math.sin(mode*(angl[l][j]))*epsilonbar;
            epsilon_der_it=-1*delta*mode*Math.sin(mode*(angl[i][t]))*epsilonbar;
            epsilon_der_ib=-1*delta*mode*Math.sin(mode*(angl[i][b]))*epsilonbar;



            lap_phi_x=((epsilon_rj*epsilon_rj)*(phi[r][j]-phi[i][j])-(epsilon_lj*epsilon_lj)*(phi[i][j]-phi[l][j]))/(delx*delx);
            lap_phi_y=((epsilon_it*epsilon_it)*(phi[i][t]-phi[i][j])-(epsilon_ib*epsilon_ib)*(phi[i][j]-phi[i][b]))/(dely*dely);

            phi_der_1=((epsilon_it*epsilon_der_it)*(phi[r][t]-phi[l][t])-(epsilon_ib*epsilon_der_ib)*(phi[r][b]-phi[l][b]))/(4*delx*dely);
            phi_der_2=((epsilon_rj*epsilon_der_rj)*(phi[r][t]-phi[r][b])-(epsilon_lj*epsilon_der_lj)*(phi[l][t]-phi[l][b]))/(4*delx*dely);
          }
          else {
            lap_phi_x=lap_phi_y=phi_der_1=phi_der_2=0;
          }

            new_phi[i][j]=phi[i][j]+(delt/tau)*(lap_phi_x + lap_phi_y + phi_der_1 - phi_der_2 + g_phi)+noise_str*phi[i][j]*(1-phi[i][j])*(-Math.random() +.5);
            term1=lheat*(new_phi[i][j]-phi[i][j]);
            term2=delt*((temp[i][t]+temp[i][b])/(dely*dely) + (temp[l][j]+temp[r][j])/(delx*delx));
            g=1+2*delt*((delx*delx)+(dely*dely))/((delx*delx)*(dely*dely));
            new_temp[i][j]= (temp[i][j] + term1 + term2)/g;

        }
    }

    for(i=0;i<nx;i++)
    {
        for(j=0;j<ny;j++)
        {
          phi[i][j]=new_phi[i][j];
          temp[i][j]=new_temp[i][j];

        }
    }

 }
}






function startStop() {
  running = !running;
  if(running){
  startButton.value = " Pause ";
}
else {
  startButton.value = " Resume ";
}
}
