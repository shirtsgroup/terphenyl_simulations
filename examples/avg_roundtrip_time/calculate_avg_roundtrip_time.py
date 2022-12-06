import heteropolymer_simulations as hs
import numpy as np
import sys
import time

def main():
    t1 = time.time()
    log_file = "../../simulations/terphenyl_mop/hexamer_remd/t_250_450/state_time_remd.png"
    log_file_obj = hs.remd_utils.REMDLogFile(log_file)
    rtts = hs.remd_utils.calculate_roundtrip_times(log_file_obj)

    print(len(rtts), "out of",log_file_obj.n_states ,"simulations complete at least 1 RT.")
    print("Average RTT: ", np.mean(rtts), "ns +/-", np.std(rtts))

    t2 = time.time()
    print("This analysis took:", round(t2 - t1), "second(s).")


if __name__ == "__main__":
    main()