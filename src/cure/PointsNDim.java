/**
 * 
 */
package cure;

public class PointsNDim {

    protected double [] coord;
    protected int index;

    protected PointsNDim(int n) {
	coord = new double [n];
    }

    protected PointsNDim(double [] x) {

	coord = new double[x.length];
	for (int i=0; i<x.length; ++i) coord[i] = x[i];
    }
    
    protected PointsNDim(double [] x, int indx) {

    	coord = new double[x.length];
    	for (int i=0; i<x.length; ++i) coord[i] = x[i];
    	index=indx;
    	
    }

    protected Object clone() {

	return new PointsNDim(coord);
    }

    protected boolean equals(PointsNDim p) {

	
	for (int i=0; i<coord.length; ++i)
	    if (coord[i] != p.coord[i])
		return false;

	return true;
    }

    protected static double sqrdist(PointsNDim x, PointsNDim y) {
	
	return PointsNDim.sqrdist(x.coord, y.coord);
    }
    
    public void copy(PointsNDim p) {
    	
    	for (int i=0; i<this.coord.length; ++i) {
    	   p.coord[i]=this.coord[i];
    	}
    	
   }

    public String toString() {
	String s = "";
	for (int i=0; i<coord.length; ++i) {
	    s = s + coord[i] + " ";
	}
	return s;
    }
    
    protected double distance(double [] a, double [] b)  {
    	
    	return Math.sqrt(sqrdist(a, b));
    	
   }
        
    protected static double sqrdist(double [] a, double [] b) {

	double dist = 0;

	for (int i=0; i<a.length; ++i) {
	    double diff = (a[i] - b[i]);
	    dist += diff*diff;
	}

	return dist;
    }   

}