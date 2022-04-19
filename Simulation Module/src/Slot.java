public class Slot {
    public double startTime;  // start time of the slot (in hours)
    public double appTime;    // appointment time of the slot, dependent on type and rule (in hours)
    public int slotType;       // type of slot (0=none, 1=elective, 2=urgent within normal working hours, 3=urgent in overtime)
    public int patientType;    // (0=none, 1=elective, 2=urgent)

    public Slot(){
        startTime = 0;
        appTime = 0;
        slotType = 0;
        patientType = 0;
    }
}
