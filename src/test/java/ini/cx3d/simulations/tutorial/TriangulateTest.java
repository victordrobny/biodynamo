package ini.cx3d.simulations.tutorial;

import ini.cx3d.BaseSimulationTest;
import ini.cx3d.spatialOrganization.SpaceNode;
import ini.cx3d.spatialOrganization.SpatialOrganizationEdge;
import ini.cx3d.spatialOrganization.SpatialOrganizationNode;

import java.util.LinkedList;
import java.util.Vector;

import static org.junit.Assert.fail;

public class TriangulateTest extends BaseSimulationTest {

    public TriangulateTest() {
        super(IntracellularDiffusionTest.class);
    }

    @Override
    public void simulate() throws Exception {

        Vector<SpatialOrganizationNode> spaceNodes = new Vector<>();
        double[] pos = {463.7047970232077, 439.8653887819098, 447.1949176631939};
        SpaceNode initialSn = new SpaceNode(pos, null);
        double[] lastPos = pos;
        long startTs = System.currentTimeMillis();
        for (int i = 0; i < 2000; i++) {
            pos = nextPos(lastPos);
            spaceNodes.add(initialSn.getNewInstance(pos, null));
            //spaceNodes.add(ecm.getPhysicalNodeInstance(randomNoise(500,3)).getSoNode());
            lastPos = pos;
        }
        long createTs = System.currentTimeMillis();
        double[] direction = {0.08741642977919392,-0.020565131563058878, -0.03049639415937795 };
        for(int i = 0; i < 100; i++) {
            for (SpatialOrganizationNode sn : spaceNodes) {
                sn.moveFrom(direction);
            }
        }
        long moveTs = System.currentTimeMillis();

        System.out.println("dt create " + (createTs - startTs));
        System.out.println("dt move   " + (moveTs - createTs));

        // validate
        int hash = 1;
        for(SpatialOrganizationNode sn : spaceNodes){
            for(SpatialOrganizationEdge se : (LinkedList<SpatialOrganizationEdge>) sn.getEdges()){
                hash = hash * 17 + ((SpaceNode)se.getOpposite(sn)).getId();
            }
        }

        // assert
        if(hash == -1570882407) {
            System.out.println("Test successful");
        } else {
            fail("wrong test result");
        }
    }

    private double[] nextPos(double[] last){
        return new double[]{ random(last[0]), random(last[1]), random(last[2]) };
    }

    private double random(double last){
        return (last * 313) % 1009;
    }
}
