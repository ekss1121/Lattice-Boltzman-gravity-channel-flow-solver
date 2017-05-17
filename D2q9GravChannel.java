import java.io.*;

public class D2q9GravChannel{

	private double[][][] f;
	private double[][][] feq;
	private double[][] rho;
	private double[][] u;
	private double[][] v;
	private double[] w;
	private double[] cx;
	private double[] cy;
	private double[] x;
	private double[] y;
	private double[] nu_t; // turb viscosity
	private double[] k;
	private double[] omega;
	public double g_lat, tau, nu_lat, nu;
	public double u_tau;
	public double cv, ch;
	public int m, n;

	public D2q9GravChannel(int reso){
		this.n = reso;
		this.m = 5;

		f = new double[m][n][9];
		feq = new double[m][n][9];
		rho = new double[m][n];
		nu_t = new double[n];
		k = new double[n];
		omega = new double[n];
		u = new double[m][n];
		v = new double[m][n];
		w = new double[9];
		cx = new double[9];
		cy = new double[9];
		x = new double[m];
		y = new double[n];
		//
		w[0] = 4./9.;
		for (int i=1; i<5; i++) w[i] = 1./9.;
		for (int i=5; i<9; i++) w[i] = 1./36.;	
	
		cx[0]=0; cx[1]=1; cx[2]=0; cx[3]=-1; cx[4]=0;
		cx[5]=1; cx[6]=-1;cx[7]=-1;cx[8]=1;
		cy[0]=0; cy[1]=0; cy[2]=1; cy[3]=0; cy[4]=-1;
		cy[5]=1; cy[6]=1; cy[7]=-1;cy[8]=-1;
		//Initialize to zero

		for (int j=0; j<n; j++){
			for (int i=0; i<m; i++){
				u[i][j] = 0;
				v[i][j] = 0;
				rho[i][j] = 1;
				nu_t[j] = 1e-7;
				k[j] = 0.01;
				omega[j] = 0.01;
				for (int k=0; k<9; k++){
					f[i][j][k] = 0.0;
					feq[i][j][k] = 0.0;
				}
			}
		}
	}


	// Pertubation
	public void pertubation(){

		double dis_y = 0;

		for (int j=1; j<n-1; j++){
			for (int i=0; i<m; i++){
				dis_y = (n-1)/2 - Math.abs(j - (n-1)/2);
				u[i][j] = (n-1)*(n-1)/(2*nu_lat)*g_lat*dis_y/(n-1)*(1-dis_y/(n-1));
			}
		}

	}


	//collision process
	public void collision(){
		double ww;
		double usquare, uv;
		double uij, vij;
		double fy = 0.0;
		double fx = g_lat;
		for (int j=0; j<n; j++){
			for (int i=0; i<m; i++){
				//ww = 1./ tau;
				ww = 1. / (3*(nu_lat+nu_t[j])+0.5);
				uij = u[i][j]; // taking acount of gravity
				vij = v[i][j];
				usquare = uij*uij + vij*vij;
				for (int k=0; k<9; k++){
					uv = uij*cx[k]+vij*cy[k];
					feq[i][j][k] = rho[i][j] * w[k]*
						(1.0+3.0*uv+4.5*uv*uv-1.5*usquare);
					f[i][j][k] = ww * feq[i][j][k]
						+(1. - ww) * f[i][j][k] + 
						 3.0 * w[k] * rho[i][j] * (cx[k]*fx + cy[k]*fy);
				}
			}
		}
	}


	// Streaming process
	public void streaming(){
		// Right to left (f1)
		for (int j=0; j<n; j++){
			for (int i=m-1; i>0; i--){
				f[i][j][1] = f[i-1][j][1];
			}
			// left to right (f3)
			for (int i=0; i<m-1; i++){
				f[i][j][3] = f[i+1][j][3];
			}
		}
		// Top to bottom (f2)
		for (int j=n-1; j>0; j--){
			for (int i=0; i<m; i++){
				f[i][j][2] = f[i][j-1][2];
			}
			// f(5)
			for (int i=m-1; i>0; i--){
				f[i][j][5] = f[i-1][j-1][5];
			}
			// f(6)
			for (int i=0; i<m-1; i++){
				f[i][j][6] = f[i+1][j-1][6];
			}
		
		}
		// Bottom to top (f4)
		for (int j=0; j<n-1; j++){
			for (int i=0; i<m; i++){
				f[i][j][4] = f[i][j+1][4];
			}
			// f(8)
			for (int i=m-1; i>0; i--){
				f[i][j][8] = f[i-1][j+1][8];
			}
			// f(7)
			for (int i=0; i<m-1; i++){
				f[i][j][7] = f[i+1][j+1][7];
			}
		}
	}

	// B.C.
	public void boundaryConds(){
		// Periodic BC on the west wall 
		for (int j=0; j<n; j++){
			f[0][j][1] = f[m-1][j][1];
			f[0][j][5] = f[m-1][j][5];
			f[0][j][8] = f[m-1][j][8];
		}
		// Periodic BC  on east wall
		for (int j=0; j<n; j++){
			f[m-1][j][3] = f[0][j][3];
			f[m-1][j][6] = f[0][j][6];
			f[m-1][j][7] = f[0][j][7];
		}
		// Bounce back on south wall
		for (int i=0; i<m; i++){
			f[i][0][2] = f[i][0][4];
			f[i][0][5] = f[i][0][7];
			f[i][0][6] = f[i][0][8];
		}
		// Bounce back on  north wall
		double rhon;
		for (int i=0; i<m; i++){
			f[i][n-1][4] = f[i][n-1][2];
			f[i][n-1][8] = f[i][n-1][6];
			f[i][n-1][7] = f[i][n-1][5];
		}
	}

	public void undatePhyFields(){
		double ssum;
		for (int j=0; j<n; j++){
			for (int i=0; i<m; i++){
				ssum = 0;
				for (int k=0; k<9; k++){
					ssum += f[i][j][k];	
				}
				rho[i][j] = ssum;
			}
		}

		//cout << rho[(M-1)*N+N-1] << endl;
		double usum;
		double vsum;		
		for (int j=1; j<n-1; j++){
			for (int i=0; i<m; i++){
				usum = 0.0;
				vsum = 0.0;
				for (int k=0; k<9; k++){
					usum += f[i][j][k]*cx[k];
					vsum += f[i][j][k]*cy[k];
				}
				u[i][j] = usum/rho[i][j];
				v[i][j] = vsum/rho[i][j];
			}
		}
	}

	public void screenInfo(int loop){
		if (loop % 5000 == 0){
			System.out.println(loop + " " + cv*u[2][(n-1)/2] + " " + ch/cv * loop);
			//System.out.println(loop + " " + nu_t[(n-1)/2] + " " + ch/cv * loop);
		}
	}

	public void init(double re_tau){
		double h = 1; // Height of the channel
		ch = h/(double)(n-1);
		double g = 10.0; // gravity
		double cs = 347.2;
		double cs_lat = 1.0 / Math.sqrt(3.0);
		cv = cs / cs_lat;
		//tau = 0.506;
		u_tau = Math.sqrt(g*h/2.0);
		nu = u_tau * h/2.0 / re_tau;
		nu_lat = nu / (cv * ch);
		tau = 3 * nu_lat + 0.5;
		double u_tau_lat = u_tau / cv;
		
		

		double umax = -g/2.0/nu*0.25 + g/2.0/nu*0.5;
		//System.out.println("nu: "+nu +" nu_lat: "+nu_lat);
		//System.out.println("Tau: "+ tau);
		//System.out.println("u_tau_lat: "+ u_tau_lat);
		g_lat = 2 * u_tau_lat * u_tau_lat / (double)(n-1);
		//System.out.println("g_lat: " + g_lat);
		for (int i=0; i<m; i++){
			x[i] = i*ch;
		}
		for (int j=0; j<n; j++){
			y[j] = j*ch;
		}
		
	}

	// Turbulence model Algebraic
	public void mixingLength(){

		double coef_k = 0.41;
		double coef_e = 0.1;
		double dh = 2 * n; // 4A/P = 4 * h * N / 2(h+N) when N -> inf
		double sij, l_mix;

		nu_t[0] = 0;
		nu_t[n-1] = 0;
		for (int j=1; j<n-1; j++){
			sij = Math.abs(0.5*(u[(m-1)/2][j+1] - u[(m-1)/2][j-1]));

			if (j < (n-1)/2){
				l_mix = Math.min(j, coef_e*dh);
			}else{
				l_mix = Math.min(n-1-j, coef_e*dh);
			}

			nu_t[j] = Math.pow(l_mix*coef_k,2.0)*sij;
		}
	}

	// Turbulence model K-Omega
	public void kOmega() {

		double sigma_k, sigma_w, bstar, r, beta;
		double k_one, k_two, nu_t_phalf, nu_t_mhalf;
		double w_one, w_two;

		sigma_k = 0.5;
		sigma_w = 0.5;
		bstar = 0.09;
		r = 5.0/9.0;
		beta = 0.75;

		nu_t[0] = 0;
		nu_t[n-1] = 0;
		k[0] = 0; k[n-1] = 0;
		omega[0] = 60*nu_lat/(beta*0.25);
		omega[n-1] = omega[0];


		for (int j=1; j<n-1; j++){
					
			nu_t_phalf = 0.5*(nu_t[j+1] + nu_t[j]);
			nu_t_mhalf = 0.5*(nu_t[j-1] + nu_t[j]);

			k_one = nu_t[j]*Math.pow(0.5*(u[(m-1)/2][j+1] - u[(m-1)/2][j-1]), 2);
			k_two = (sigma_k*(nu_lat+nu_t_phalf)*(k[j+1]-k[j]) - 
				sigma_k * (nu_lat + nu_t_mhalf)*(k[j]-k[j-1]));

			w_one = r*Math.pow(0.5*(u[(m-1)/2][j+1] - u[(m-1)/2][j-1]), 2);
			w_two = (sigma_w*(nu_lat+nu_t_phalf)*(omega[j+1]-omega[j]) - 
				sigma_w * (nu_lat + nu_t_mhalf)*(omega[j]-omega[j-1]));

			k[j] = k[j] + k_one + k_two - bstar*omega[j]*k[j];
			omega[j] = omega[j] + w_one + w_two - beta* omega[j] * omega[j];

			nu_t[j] = k[j] / omega[j];

		}
	}

	public void writeData() throws IOException{

		FileWriter vmag = null;
		FileWriter vmid = null;
		FileWriter upyp = null;
		 
		double velMag;
		double nu_total;
		try{
			vmag = new FileWriter("vmag.dat");
			vmid = new FileWriter("umid-RE20.dat");
			upyp = new FileWriter("uplus-vs-yplus.dat");
			vmag.write("# variables = x[m] y [m]  vel_magnitude [m/s] \n");
			vmid.write("# variables = y [m]  x-velocity [m/s] \n");
			upyp.write("# variables = y+ u+ \n");
			for (int j=0; j<n; j++){
				for (int i=0; i<m; i++){
					velMag = cv * Math.sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
					vmag.write(x[i] + " " +y[j] + " " + u[i][j] + "\n");
					if (i == 2){
						vmid.write(y[j] + " " + cv * u[i][j] + "\n");
						nu_total = (nu_lat + 0) * cv * ch;
						upyp.write(y[j] * u_tau /nu_total + " " + u[i][j] * cv / u_tau+"\n");
					}
				}
				vmag.write("\n");
			}
		} finally{
			if (vmag != null){
				vmag.close();
			}
			if (vmid != null){
				vmid.close();
			}
			if (upyp != null){
				upyp.close();
			}
		}
		
		 
		
	}
	
}