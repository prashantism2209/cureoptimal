package cure;
import java.util.Comparator;

/**
 * Defines a class CompareCluster which helps the MinHeap (Implemented using Priority Queue) 
 * to compare two clusters and store accordingly in the heap.
 * 
 * The 2 clusters are compared based on the distance from their closest Cluster. The cluster pair which has the lowest such distance
 * is stored at the root of the min heap.
 
 */
public class CompareCluster implements Comparator{

	public int compare(Object CLUSTER1, Object CLUSTER2) {
		Cluster cluster1 = (Cluster)CLUSTER1;
		Cluster cluster2 = (Cluster)CLUSTER2;
		if(cluster1.distanceFromClosest < cluster2.distanceFromClosest) {
			return -1;
		}
		else if(cluster1.distanceFromClosest == cluster2.distanceFromClosest) {
			return 0;
		}
		else return 1;
	}
}
