import java.io.*;

public class GravSims{


	public static void main(String[] args){
		int iters = Integer.parseInt(args[0]);
		double re_tau = Double.parseDouble(args[1]);
		int reso = 151;
		D2q9GravChannel mySims1 = new D2q9GravChannel(reso);
		//Calculate lattice parameters
		// Generate computing grids;
		mySims1.init(re_tau);
		//mySims1.pertubation();
		// Main loop

		for (int loop=1; loop<=iters; loop++){
			//mySims1.kOmega();
			//mySims1.mixingLength();
			mySims1.collision();
			mySims1.streaming();
			mySims1.boundaryConds();
			mySims1.undatePhyFields();
			mySims1.screenInfo(loop);
				
		}

		try{
			mySims1.writeData();
		} catch (IOException e){
			System.out.println("General I/O Exception: " + e.getMessage());
		}

		
	}
}