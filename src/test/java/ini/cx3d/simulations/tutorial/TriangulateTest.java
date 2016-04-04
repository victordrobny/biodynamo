package ini.cx3d.simulations.tutorial;

import ini.cx3d.BaseSimulationTest;
import ini.cx3d.spatialOrganization.factory.SpaceNodeFactory;
import ini.cx3d.spatialOrganization.interfaces.SpaceNode;
import ini.cx3d.spatialOrganization.SpatialOrganizationEdge;
import ini.cx3d.spatialOrganization.SpatialOrganizationNode;
import ini.cx3d.swig.spatialOrganization.ListT_Edge;

import java.util.LinkedList;
import java.util.Vector;

import static org.junit.Assert.fail;

public class TriangulateTest extends BaseSimulationTest {

    public TriangulateTest() {
        super(IntracellularDiffusionTest.class);
    }

    @Override
    public void simulate() throws Exception {

        Vector<ini.cx3d.spatialOrganization.interfaces.SpaceNode> spaceNodes = new Vector<>();
        double[] pos = {463.7047970232077, 439.8653887819098, 447.1949176631939};
        SpaceNode initialSn = new SpaceNodeFactory().create(pos, null);
        double[] lastPos = pos;
        long startTs = System.currentTimeMillis();
        for (int i = 0; i < 2000; i++) {
            pos = nextPos(lastPos);
            spaceNodes.add((SpaceNode) initialSn.getNewInstance(pos, null));
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
            for(SpatialOrganizationEdge se : (ListT_Edge) sn.getEdges()){
                hash = hash * 17 + ((SpaceNode)se.getOpposite((SpaceNode) sn)).getId();
            }
        }

        // assert
        if(hash == -1669558583) {
            System.out.println("Test successful");
        } else {
            System.out.println(hash);
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
