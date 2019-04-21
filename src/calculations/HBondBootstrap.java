package calculations;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;

/**
 * Used to estimate confidence intervals of HB dynamic and diffusive parameters
 * @author trevorfisher
 *
 */
public class HBondBootstrap 
{
	private double confInterval;
	private double kbStart;
	private double krStart;
	private double taoStart;
	private double c0;
	public HBondBootstrap(double kb, double kr, double tao, double c0, double interval) {
		//Init. starting values, conf. interval
		this.kbStart = kb;
		this.krStart = kr;
		this.taoStart = tao;
		this.c0 = c0;
		this.confInterval = interval;
	}
	
	public void runAnalysis() {
		Scanner sc = new Scanner(System.in);
		System.out.println("Provide output file for HB bootstrapping.");
		String fileName = sc.nextLine();
		runAnalysis(fileName);
	}
	
	public void runAnalysis(String fileName) {
		//Read in file containing experimental c(t) and conf. intervals
		ArrayList<Double> times = new ArrayList<Double>();
		ArrayList<Double> cTExp = new ArrayList<Double>();
		ArrayList<Double> cTConf = new ArrayList<Double>();
		parseOutput(fileName, times, cTExp, cTConf);
		//System.out.println(times.size() + " " + cTExp.size() + " " + cTConf.size());
		//Store values of model that can fit more than the conf. interval of data.
		ArrayList<Double> kb = new ArrayList<Double>();
		ArrayList<Double> kr = new ArrayList<Double>();
		ArrayList<Double> tao = new ArrayList<Double>();
		
		double tails = 0.65;
		double strikePerc = confInterval;
		double deviationkb = .2;
		double deviationkr = .2;
		double deviationtao = .4;
		boolean converged = false;
		while(!converged) {
			while(kb.size() < 400) {
				//Generate model values
				//Allow value to be generated that lies a deviation above or below
				double newKB = -1, newKR = -1, newTAO = -1;
				if(deviationkb < 1.0) {
					newKB = kbStart+(2*Math.random()-1.0)*deviationkb*kbStart;
				} else {
					newKB = deviationkb*kbStart+(2*Math.random()-1.0)*deviationkb*kbStart;
				}
				if(deviationkr < 1.0) {
					newKR = krStart+(2*Math.random()-1.0)*deviationkr*krStart;
				} else {
					newKR = deviationkr*krStart+(2*Math.random()-1.0)*deviationkr*krStart;
				}
				if(deviationtao < 1.0) {
					newTAO = taoStart+(2*Math.random()-1.0)*deviationtao*taoStart;
				} else {
					newTAO = deviationtao*taoStart+(2*Math.random()-1.0)*deviationtao*taoStart;
				}
				//System.out.println(newKB + " " + newKR + " " + newTAO);
				int strikes = 0;
				for(int i = 0; i < cTExp.size(); i++) {
					//Compare model cT with exp cT, see if within conf. interval
					double cTModel = Stehfest_ct(times.get(i), newKB, newKR, newTAO, c0);
					//Add a strike if it doesn't fit within conf. intervals at this time point
					if(cTModel > cTExp.get(i)+cTConf.get(i) || cTModel < cTExp.get(i)-cTConf.get(i)) {
						strikes += 1;
						//System.out.println("Strike");
					}
					//Once enough times have been violated, end
					if(strikes > cTExp.size()*(1-strikePerc)) {
						break;
					}
				}
				if(strikes <= cTExp.size()*(1-strikePerc)) {
					//Store model values that fit the data
					kb.add(newKB);
					kr.add(newKR);
					tao.add(newTAO);
					//System.out.println(newKB + " " + newKR + " " + newTAO);
					if(kb.size()%20 == 0) {
						System.out.println(kb.size());
					}
				} else {
					//System.out.println("Failed combo");
				}
			}
			//Check to see if distribution is wide enough
			Collections.sort(kb);
			Collections.sort(kr);
			Collections.sort(tao);
			
			double range1 = kb.get((int)(kb.size()*(1-(1-confInterval)/2)))-kb.get((int)(kb.size()*(1-confInterval)/2));
			double range2 = kr.get((int)(kr.size()*(1-(1-confInterval)/2)))-kr.get((int)(kr.size()*(1-confInterval)/2));
			double range3 = tao.get((int)(tao.size()*(1-(1-confInterval)/2)))-tao.get((int)(tao.size()*(1-confInterval)/2));
			double sumkb = 0; double sumkr = 0; double sumtao = 0;
			for(int i = (int)(kb.size()*(1-confInterval)/2); i < (int)(kb.size()*(1-(1-confInterval)/2)); i++) {
				sumkb += kb.get(i);
				sumkr += kr.get(i);
				sumtao += tao.get(i);
			}
			sumkb = sumkb / kb.size() / confInterval;
			sumkr = sumkr / kr.size() / confInterval;
			sumtao = sumtao / tao.size() / confInterval;
			
			if(range1 < 2*kbStart*deviationkb*tails && range2 < 2*krStart*deviationkr*tails && range3 < 2*taoStart*deviationtao*tails){
				converged = true;
			} else {
				double bench1 = range1/2/kbStart/deviationkb;
				double bench2 = range2/2/krStart/deviationkr;
				double bench3 = range3/2/taoStart/deviationtao;
				System.out.println("---");
				System.out.println("Ranges: " + bench1 + " " + bench2 + " " + bench3);
				if(bench1-tails > 0) {
					deviationkb = deviationkb + (bench1-tails+0.05)*3*deviationkb;
				} else if(tails-bench1 > 0.2) {
					deviationkb = 0.7*deviationkb;
				}
				if(bench2-tails > 0) {
					deviationkr = deviationkr + (bench2-tails+0.05)*3*deviationkr;
				} else if(tails-bench2 > 0.2) {
					deviationkr = 0.7*deviationkr;
				}
				if(bench3-tails > 0) {
					deviationtao = deviationtao + (bench3-tails+0.05)*3*deviationtao;
				} else if(tails-bench3 > 0.2) {
					deviationtao = 0.7*deviationtao;
				}
				System.out.println("Adjusted deviations: " + deviationkb + " " + deviationkr + " " + deviationtao);
				kbStart = sumkb;
				krStart = sumkr;
				taoStart = sumtao;
				System.out.println("New averages: " + kbStart + " " + krStart + " " + taoStart);
				System.out.println("---");
				kb.clear();
				kr.clear();
				tao.clear();
			}
		}
		//Print averages and deviations
		double rangekb = kb.get((int)(kb.size()*0.95))-kb.get((int)(kb.size()*0.05));
		double rangekr = kr.get((int)(kr.size()*0.95))-kr.get((int)(kr.size()*0.05));
		double rangetao = tao.get((int)(tao.size()*0.95))-tao.get((int)(tao.size()*0.05));
		double sumkb = 0; double sumkr = 0; double sumtao = 0;
		for(int i = (int)(kb.size()*(1-confInterval)/2); i < (int)(kb.size()*(1-(1-confInterval)/2)); i++) {
			sumkb += kb.get(i);
			sumkr += kr.get(i);
			sumtao += tao.get(i);
		}
		sumkb = sumkb / kb.size() / confInterval;
		sumkr = sumkr / kr.size() / confInterval;
		sumtao = sumtao / tao.size() / confInterval;
		double mediankb = kb.get(kb.size()/2);
		double mediankr = kr.get(kr.size()/2);
		double mediantao = tao.get(tao.size()/2);
		
		System.out.println("Value: average/median +- " + confInterval + "% conf.");
		System.out.println("kb: " + sumkb + "/" + mediankb + " +- " + rangekb/2);
		System.out.println("kr: " + sumkr + "/" + mediankr + " +- " + rangekr/2);
		System.out.println("tao: " + sumtao + "/" + mediantao + " +- " + rangetao/2);
		
		System.out.println("kb min: " + kb.get((int)(kb.size()*0.05)) + " - kb max: " + kb.get((int)(kb.size()*0.95)));
		System.out.println("kr min: " + kr.get((int)(kr.size()*0.05)) + " - kr max: " + kr.get((int)(kr.size()*0.95)));
		System.out.println("tao min: " + tao.get((int)(tao.size()*0.05)) + " - tao max: " + tao.get((int)(tao.size()*0.95)));
	}
	
	/**
	 * A numerical inverse Laplace transform for diffusion-paired HB dynamics
	 * @param time
	 * @param kb
	 * @param kr
	 * @param tao
	 * @param c0
	 * @return
	 */
	public static double Stehfest_ct(double time, Double kb, Double kr, Double tao, Double c0) {
		if(time == 0) {
			return c0;
		}
		int n = 16;
		int hn = n/2;
		double sum = 0.0;
		for(int i = 1; i <= n; i++) {
			//Compute value for Vi
			double vsum = 0.0;
			for(int k = (int)((i+1)/2); k <= Math.min(i, hn); k++) {
				vsum = vsum + (Math.pow(k, hn) * factorial(2*k)) / 
						(factorial(hn - k) * factorial(k) * factorial(k-1) * factorial(i-k) * factorial(2*k-i));
			}
			vsum = vsum * Math.pow(-1, (hn+i));
			//Multiply Vi by Fs_ct(i) term
			sum = sum + vsum * Fs_ct(Math.log(2)*i/time, kb, kr, tao, c0);
		}
		sum = sum * Math.log(2) / time;
		return sum;
	}
	
	public static double Fs_ct(double s, double kb, double kr, double tao, double c0) {
		double val = -kb * c0 / (s * (s + kb + kr * s * f_1(s, tao))) + c0 / s;
		return val;
	}
	
	public static double f_1(double s, double tao) {
		double val = 3 * tao * (1-Math.sqrt(s*tao)*Math.atan(1/(Math.sqrt(s*tao))));
		return val;
	}
	
	public static long factorial(int n) {
		long result = 1;
		for(int i = n; i > 1; i--) {
			result = result * i;
		}
		return result;
	}
	
	private static void parseOutput(String fileName, ArrayList<Double> times, ArrayList<Double> cTExp, ArrayList<Double> cTConf) {
		try {
				File file = new File(fileName);
				FileReader fr = new FileReader(file);
				BufferedReader br = new BufferedReader(fr);
				//Read first line in to have same offset
				//br.readLine();
				String nextLine;
				nextLine = br.readLine();
				while(nextLine != null && !nextLine.isEmpty()) {
					String parts[] = nextLine.split("\\s+");
					times.add(Double.valueOf(parts[0]));
					cTExp.add(Double.valueOf(parts[1]));
					cTConf.add(Double.valueOf(parts[2]));
					nextLine = br.readLine();
				}
			} catch (IOException ex) {
				ex.getMessage();
			}
	}
}
