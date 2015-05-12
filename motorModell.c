/*-----------------------------------------------------------------

Author: Max Schachtschabel
Version: 0.4

Differential equations for a simple motor model.
-------------------------------------------------------------------
*/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.141592653589793238462643383279502884

struct param {
	
	double Vd;	//Versorgungsspannung
	double R;	//Phasenwiderstand
	double L;	//Phaseninduktivität
	double kw;	//Spannungs-Konstante
	double kt;	//Drehmoment-Konstante
	double J;	//Trägheitsmoment(Motor)
	double B;	//Reibungskoeffizient
	double p;	//Anzahl der Pole
	double r; 	//parms->rius des Motors
	double m;	//Masse des Fahrzeugs
	double cwA;	//Windwiderstandswert
	double crr;	//Reibungskoeffizient	
	double rho;	
};


struct input {
	double T_load;
};


struct output {
	double t;
	double Va;
	double Vb;
	double Vc;
	double bEMFa;
	double bEMFb;
	double bEMFc;
	double ia;
	double ib;
	double ic;
	double phi;
	double v;
	double Te;
	double s;
	double a;
	double P;
	int sig;
};

struct vector {
	double t;
	double ia;
	double ib;
	double ic;
	double s;
	double w;
	double a;
	double P;
};

FILE *fd;

void print_out (struct output *out, double t) {
    //fprintf(fd,"%f \t",out->t);
    fprintf(fd,"%f \t",out->s);	
    //printf("%f \t %f \t %f \t",out->Va,out->Vb,out->Vc);
    //printf("%f \t %f \t %f \t",out->ia,out->ib,out->ic);
    fprintf(fd,"%f\n",out->v);
    //fprintf(fd,"%f\n",out->P);
    //printf("%f\t",out->phi);
    //printf("%f \t %f \t %f \t",out->bEMFa,out->bEMFb,out->bEMFc);
    //printf("%f\t",out->Te);
    //printf("%f \n",out->P);
}


void scale (struct vector *v,double h) {
	v->t *= h;
	v->ia *= h;
	v->ib *= h;
	v->ic *= h;
	v->s *= h;
	v->w *= h;
	v->P *= h;
}


void copy_Vec (struct vector *from, struct vector *to) {
	to->t = from->t;	
	to->ia = from->ia;
	to->ib = from->ib;
	to->ic = from->ic;
	to->s = from->s;
	to->w = from->w;
	to->P = from->P;
}


void add_Vec (struct vector *to, struct vector *from) {
	to->t += from->t;	
	to->ia += from->ia;
	to->ib += from->ib;
	to->ic += from->ic;
	to->s += from->s;
	to->w += from->w;
	to->P += from->P;	
}


void calc_diff (struct vector *v, struct input *in, struct output *out, struct param *parms, struct vector *deriv) {

	double mphi,phih,rest,v_ms,phiA,phiB,phiC;

	if (out->phi >= 360) { 
		out->phi -= 360; 
	}

	phih = out->phi*(parms->p/2.0);
	
	if (phih<0) phih = 0;
	rest = phih - (int)phih;
	phih = ((int)phih % 360) + rest;

	mphi = phih;
	rest = mphi - (int)mphi;
	mphi = ((int)mphi % 60) + rest; 

	if (phih <= 60) {
		out->Va = parms->Vd/2; out->Vb = 0; out->Vc= -parms->Vd/2;
		phiA = -1+mphi/30; phiB = 1; phiC = -1;} 
	else if (phih <= 120) {
		out->Va = parms->Vd/2; out->Vb = -parms->Vd/2; out->Vc = 0;
		phiA = 1; phiB = 1-mphi/30; phiC = -1;}
	else if (phih <= 180) {
		out->Va = 0; out->Vb = -parms->Vd/2; out->Vc = parms->Vd/2;
		phiA = 1; phiB = -1; phiC = -1+mphi/30;}
	else if (phih <= 240) {
		out->Va = -parms->Vd/2; out->Vb = 0; out->Vc = parms->Vd/2;
		phiA = 1-mphi/30; phiB = -1; phiC = 1;}
	else if (phih <= 300) {
		out->Va = -parms->Vd/2; out->Vb = parms->Vd/2; out->Vc = 0;
		phiA = -1; phiB = -1+mphi/30; phiC = 1;}
	else {
		out->Va = 0; out->Vb = parms->Vd/2; out->Vc = -parms->Vd/2;
		phiA = -1; phiB = 1; phiC = 1-mphi/30;
	}

	if (out->sig == 0) { out->Va = 0; out->Vb = 0; out->Vc = 0; };

	v_ms = v->w * (parms->r/(2*PI));
	
	double e_a = parms->kw*phiA*v->w; 
	double e_b = parms->kw*phiB*v->w; 
	double e_c = parms->kw*phiC*v->w;

	double T_a = parms->kt*phiA*v->ia; 
	double T_b = parms->kt*phiB*v->ib;
	double T_c = parms->kt*phiC*v->ic;

	double M_mot = T_a + T_b + T_c - in->T_load - parms->B*v->w;
	double F_drag = 0.5*parms->rho*parms->cwA*v_ms*v_ms + parms->m*9.81*parms->crr + parms->m*9.81*0;

	/* Derivate w.R.t. Time 
	deriv->ia = (out->Va - e_a - parms->R*v->ia) / parms->L;
    	deriv->ib = (out->Vb - e_b - parms->R*v->ib) / parms->L;
    	deriv->ic = (out->Vc - e_c - parms->R*v->ic) / parms->L;
	deriv->w = (M_mot - F_drag*parms->r) / (parms->J + parms->m * parms->r * parms->r);
	deriv->t = 1;
	deriv->s = v_ms;
	deriv->P = (M_mot / parms->r) * v_ms; 	
	*/

	/* Derivate w.R.t. Distance */
	deriv->ia = (out->Va - e_a - parms->R*v->ia) / (parms->L * v_ms);
	deriv->ib = (out->Vb - e_b - parms->R*v->ib) / (parms->L * v_ms);
    	deriv->ic = (out->Vc - e_c - parms->R*v->ic) / (parms->L * v_ms);
	deriv->w = (M_mot - F_drag*parms->r) / ((parms->J + parms->m * parms->r * parms->r) * v_ms);
	deriv->t = 1/v_ms;
	deriv->s = 1;
	deriv->P = (M_mot / parms->r) < 0 ? 0 : M_mot / parms->r;

	out->a = deriv->w;
	out->bEMFa = e_a;
	out->bEMFb = e_b;
	out->bEMFc = e_c;

	//printf("%.30f,%.30f,%.30f,%.30f,%.30f\n",v_ms,deriv->P,M_mot,F_drag,phih);
}


void step_diff (struct vector *v, struct input *in, struct output *out,struct param *parms, double h) {

    	struct vector k[4],t[4];
		
	calc_diff(v,in,out,parms,&k[0]);
	scale(&k[0],h);
	copy_Vec(&k[0],&t[0]);
	scale(&t[0],0.5);
	add_Vec(&t[0],v);

	calc_diff(&t[0],in,out,parms,&k[1]);
	scale(&k[1],h);
	copy_Vec(&k[1],&t[1]);
	scale(&t[1],0.5);
	add_Vec(&t[1],v);

	calc_diff(&t[1],in,out,parms,&k[2]);
	scale(&k[2],h);
	copy_Vec(&k[2],&t[2]);
	add_Vec(&t[2],v);

	calc_diff(&t[2],in,out,parms,&k[3]);
	scale(&k[3],h);

	scale(&k[1],2);
	scale(&k[2],2);
	add_Vec(&k[0],&k[1]);
	add_Vec(&k[0],&k[2]);
	add_Vec(&k[0],&k[3]);
	scale(&k[0],1.0/6.0);
	add_Vec(v,&k[0]);
} 





void motor_eq (double t, struct input *in, struct param *parms, double h) {

	double i,w_last,t_last,phih,rest,P,*PP,mean,change,P_last;
	int j,k,Plen,interv,runs=0;
	
	struct output out;
	struct vector v_next;
	
	interv = 10;
	Plen = (int)(t*(1/h));
	out.sig = 1;
	change = 100;	

	v_next.ia = 0;
	v_next.ib = 0;
	v_next.ic = 0;
    	v_next.w = 2 / (parms->r/(2*PI));
	out.P = 0;

	for (i=0;i<t;i+=h) {
		
		w_last = v_next.w;
		t_last = v_next.t;
		P_last = out.P;
/*
		if (i > change) {
			out.sig = (out.sig+1)%2;
			change += 100;
		}
*/

        	step_diff(&v_next,in,&out,parms,h);

		out.phi += h*(v_next.w + w_last)/2.0;

		out.ia = v_next.ia;
        	out.ib = v_next.ib;
        	out.ic = v_next.ic;
		out.v = v_next.w * (parms->r/(2*PI));
		out.Te = out.a * parms->J;
		out.s = v_next.s;
		out.P = v_next.P;

		if ((runs % 1000) == 0) {
			print_out(&out,t);
		}
		runs++;
	}
}


void setparms (struct param *parms, float var_parms[]) {

	parms->Vd = 40.0;
	parms->p = 42.0; 
	parms->r = 0.23;
	parms->m = 100.0;
	parms->cwA = 0.2;
	parms->crr = 0.01;
	parms->rho = 1.2;

	parms->R = 0.089;
	parms->L = 0.000392;
	parms->kw = 0.045;	
	parms->kt = 0.75;
	parms->J = 0.102;
	parms->B = 0.01;
/*
	parms->R = var_parms[0];
	parms->L = var_parms[1];
	parms->kw = var_parms[2];	
	parms->kt = var_parms[3];
	parms->J = var_parms[4];
	parms->B = var_parms[5];
*/
}


void main (int argc, char* argv[]) {

	struct param parms;
	struct input in;
	int t,iter,i=0;
	int pos = 0;
	double h = 0.00025;
	float var_parms[6];
	float parm;

    	t = atoi(argv[1]);	
	FILE *parmfile = fopen(argv[2],"r");


	while(fscanf(parmfile,"%f\n",&parm) != EOF) {	
		var_parms[i] = parm;
		printf("%f\n",var_parms[i]);
		i++;
	}

	setparms(&parms,var_parms);

	in.T_load = 0;
	fd = fopen("model_data","w");
	motor_eq(t,&in,&parms,h);
	return;	

}



