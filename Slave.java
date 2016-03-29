/**
 * Created by john on 3/29/2016.
 */
public class Slave extends Thread {
    public void run() {
        System.out.println("Slave " + this.getId());
    }
    public static void main(String[] args) {
        Slave obj = new Slave();
        obj.start();
    }
}
