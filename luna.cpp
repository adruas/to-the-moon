/*Programa que simula la interacción de una nave enviada desde la Tierra con la Luna*/

//BIBLIOTECAS
#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;

//CONSTANTES
const int N=4; //Numero de ecuaciones diferenciales
const double omega=2.6617e-6, delta=7.018782618e-12, mu=0.01230246418; //Constantes

//PROTOTIPOS DE FUNCIONES
double radio(double x1, double x2, double x3, double x4); //Función radial
double ang(double y1, double y2, double y3, double y4); //Función angular
double momr(double z1, double z2, double z3, double z4, double tiempo1); //Función momento radial
double momang(double t1, double t2, double t3, double t4, double tiempo2); //Función momento angular
void ca1 (double (&k1)[N], double R, double PHI, double PR, double PPHI, double H, double tiempo); //Función parámetro k1
void ca2 (double (&k2)[N], double k1[], double R, double PHI, double PR, double PPHI, double H, double tiempo); //Función parámetro k2
void ca3 (double (&k3)[N], double k2[], double R, double PHI, double PR, double PPHI, double H, double tiempo);//Función parámetro k3
void ca4 (double (&k4)[N], double k3[], double R, double PHI, double PR, double PPHI, double H, double tiempo);//Función parámetro k4
double Hamiltoniano(double perre, double erre, double pefi, double fifi, double tempo); //Funcion hamiltoniano

/****************
FUNCIÓN PRINCIPAL
*****************/
int main(void)
{
 //Declaración de variables
 double t, h, r0, phi0, p_r0, p_phi0, K1[N], K2[N], K3[N], K4[N];
 double r,phi,p_r,p_phi, x_nave, y_nave, x_luna, y_luna , Hamilton;
 int l;
 ofstream luna, nave, hamiltoniano; //Ficheros

 //Valores para h y semillas (posiciones iniciales y momentos iniciales)
 cout << "Introduzca un valor de h: ";
 cin >> h;
 r0=0.016592507804370448; phi0=0.5435; p_r0=2.9e-5; p_phi0=3.2e-15;

 /*-------ALGORITMO DE RUNGE-KUTTA-------*/
 cout << "Introduzca el número de pasos del algoritmo: ";
 cin >> l;
 //Inicializo las coordenadas:
 r=r0; phi=phi0; p_r=p_r0; p_phi=p_phi0;
 t=0.0;

 //Apertura de ficheros
 luna.open("trayectoria_luna.txt");
 nave.open("trayectoria_nave.txt");
 hamiltoniano.open("hamiltoniano.txt");
 if(!(luna.is_open()||nave.is_open()||hamiltoniano.is_open())) cout <<"\n\nError al abrir uno de los ficheros\n\n";
 luna << "#Posición x \t Posición y" << endl;
 luna << cos(omega*t) << "\t" << sin(omega*t) << endl;
 nave << "#Posición x \t Posición y" << endl;
 nave << r*cos(phi) << "\t" << r*sin(phi) << endl;
 hamiltoniano << "#Tiempo \t Hamiltoniano" << endl;
 hamiltoniano << t << "\t" << Hamiltoniano(p_r, r, p_phi, phi, t) << endl;

 //Bucle
 while(t<l)
 {
  //Creo los vectores K1, K2, K3 y K4
  ca1(K1,r,phi,p_r,p_phi,h,t); ca2(K2,K1,r,phi,p_r,p_phi,h,t);
  ca3(K3,K2,r,phi,p_r,p_phi,h,t); ca4(K4,K3,r,phi,p_r,p_phi,h,t);
  //Actualizo los valores de las coordenadas
  r=r+((1.0/6.0)*(K1[0]+2*K2[0]+2*K3[0]+K4[0]));
  phi=phi+((1.0/6.0)*(K1[1]+2*K2[1]+2*K3[1]+K4[1]));
  if(int(t)%997==0) nave << r*cos(phi) << "\t" << r*sin(phi) << endl; //Elijo 997 porque es primo
  p_r=p_r+((1.0/6.0)*(K1[2]+2*K2[2]+2*K3[2]+K4[2]));
  p_phi=p_phi+((1.0/6.0)*(K1[3]+2*K2[3]+2*K3[3]+K4[3]));
  if(int(t)%997==0) luna << cos(omega*t) << "\t" << sin(omega*t)<<endl; //Elijo 997 porque es primo
  Hamilton=Hamiltoniano(p_r, r, p_phi, phi, t);
  if(int(t)%71==0) hamiltoniano << t << "\t" << Hamilton << endl; //Elijo 5 porque es primo y necesito mas puntos
  //Avanzo un paso más
  t=t+h;
 }
 luna.close(); nave.close();
 return 0;
}



/*******************************************************/
/*FUNCIONES PARA LAS DISTINTAS ECUACIONES DIFERENCIALES*/
/*******************************************************/
//OJO: _1 es r, _2 es phi, _3 es p_r, _4 es p_phi
/*Función radial*/
double radio(double x3)
{
 return x3;
}
/*Función angular*/
double ang(double y1,double y4)
{
 return y4/(y1*y1);
}
/*Función momento radial*/
double momr(double z1, double z2, double z4, double tiempo1)
{
 double erreprima;
 erreprima=sqrt(1+(z1*z1)-(2.0*z1*cos(z2-(omega*tiempo1))));
 return (z4*z4)/(z1*z1*z1)-delta*((1.0/(z1*z1))+((mu/(erreprima*erreprima*erreprima))*(z1-cos(z2-(omega*tiempo1)))));
}
/*Función momento angular*/
double momang(double t1, double t2, double tiempo2)
{
 double errepr=sqrt(1+(t1*t1)-(2.0*t1*cos(t2-(omega*tiempo2))));
 return -delta*mu*t1*sin(t2-(omega*tiempo2))/(errepr*errepr*errepr);
}



/***********************************************/
/*FUNCIONES PARA LAS COMPONENTES DE LA MATRIZ K*/
/***********************************************/
/*Función para el cálculo de vector K1 de r, phi, p_r, p_phi*/
void ca1 (double (&k1)[N], double R, double PHI, double PR, double PPHI, double H, double tiempo)
{
 k1[0]=H*radio(PR);
 k1[1]=H*ang(R,PPHI);
 k1[2]=H*momr(R,PHI,PPHI,tiempo);
 k1[3]=H*momang(R,PHI,tiempo);
 return;
}

/*Función para el cálculo de vector K2 de r, phi, p_r, p_phi*/
void ca2 (double (&k2)[N], double k1[], double R, double PHI, double PR, double PPHI, double H, double tiempo)
{
 k2[0]=H*radio(PR+0.5*k1[2]);
 k2[1]=H*ang(R+0.5*k1[0],PPHI+0.5*k1[3]);
 k2[2]=H*momr(R+0.5*k1[0],PHI+0.5*k1[1],PPHI+0.5*k1[3],tiempo+0.5*H);
 k2[3]=H*momang(R+0.5*k1[0],PHI+0.5*k1[1],tiempo+0.5*H);
 return;
}

/*Función para el cálculo de vector K3 de r, phi, p_r, p_phi*/
void ca3 (double (&k3)[N], double k2[], double R, double PHI, double PR, double PPHI, double H, double tiempo)
{
 k3[0]=H*radio(PR+0.5*k2[2]);
 k3[1]=H*ang(R+0.5*k2[0],PPHI+0.5*k2[3]);
 k3[2]=H*momr(R+0.5*k2[0],PHI+0.5*k2[1],PPHI+0.5*k2[3],tiempo+0.5*H);
 k3[3]=H*momang(R+0.5*k2[0],PHI+0.5*k2[1],tiempo+0.5*H);
 return;
}

/*Función para el cálculo de vector K4 de r, phi, p_r, p_phi*/
void ca4 (double (&k4)[N], double k3[], double R, double PHI, double PR, double PPHI, double H, double tiempo)
{
 k4[0]=H*radio(PR+k3[2]);
 k4[1]=H*ang(R+k3[0],PPHI+k3[3]);
 k4[2]=H*momr(R+k3[0],PHI+k3[1],PPHI+k3[3],tiempo+H);
 k4[3]=H*momang(R+k3[0],PHI+k3[1],tiempo+H);
 return;
}

/*Función para el cálculo del hamiltoniano*/
double Hamiltoniano(double perre, double erre, double pefi, double fifi, double tempo)
{
  double Ham, rprima;
  rprima=sqrt(1+(erre*erre)-(2.0*erre*cos(fifi-(omega*tempo))));
  Ham=(perre/2.0)-(1.0*delta/erre)-(1.0*delta*mu/rprima)+(pefi/(2.0*erre*erre))-(omega*pefi);
  return Ham;
}
