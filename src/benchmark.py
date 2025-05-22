import os
import subprocess
import time


class NetworkConfiguration:
    
    def __init__(self, id, bandwidth, delay, jitter=None, packet_loss=None,
                 burst_size=None, quantum=None, r2q=None):
        """
        Create a network configuration.

        Parameters:
        bandwidth (int): The maximum bandwidth limit in megabits per second (Mbps).
        delay (int): The network delay to apply, in milliseconds.
        jitter (int, optional): Variability in network delay, in milliseconds.
        loss (float, optional): Packet loss rate as a percentage (e.g., 0.1 for 0.1%).
        burst (int, optional): Burst size in kilobytes, allowing short bursts over the rate limit.
        quantum (int, optional): Quantum size in bytes.
        r2q (int, optional): Rate to quantum ratio.
        """
        self.id = id
        self.bandwidth = bandwidth
        self.delay = delay
        self.jitter = jitter
        self.packet_loss = packet_loss
        self.burst_size = burst_size
        self.quantum = quantum
        self.r2q = r2q

    def __str__(self):
        return (
        f"Network configuration:\n"
        f"    id............{self.id if self.id is not None else 'N/A'}\n"
        f"    bandwidth.....{self.bandwidth if self.bandwidth is not None else 'N/A':<6} (Mbps)\n"
        f"    delay.........{self.delay if self.delay is not None else 'N/A':<6} (ms)\n"
        f"    jitter........{self.jitter if self.jitter is not None else 'N/A':<6} (ms)\n"
        f"    packet_loss...{self.packet_loss if self.packet_loss is not None else 'N/A':<6} (%)\n"
        f"    burst_size....{self.burst_size if self.burst_size is not None else 'N/A':<6} (KB)\n"
        f"    quantum.......{self.quantum if self.quantum is not None else 'N/A':<6} (bytes)\n"
        f"    r2q...........{self.r2q if self.r2q is not None else 'N/A':<6}\n"
    )


default_network_configurations = {
    "LAN (high)"                      : NetworkConfiguration("LAN10000",    10000, 0.1,   0.01,  None,  None,  9000,  25),
    "LAN (medium)"                    : NetworkConfiguration("LAN1000",     1000,  0.3,   0.02,  None,  None,  3000,  15),
    "LAN (low)"                       : NetworkConfiguration("LAN500",      500,   1,     0.2,   None,  None,  1500,  10),
    "Corporate WAN (high resilience)" : NetworkConfiguration("crpWAN500",   500,   50,    5,     0.1,   1000,  1500,  15),
    "Regional WAN (high)"             : NetworkConfiguration("regWAN500",   500,   10,    2,     0.1,   1500,  1500,  25),
    "Regional WAN (medium)"           : NetworkConfiguration("regWAN250",   250,   15,    5,     0.1,   1000,  1500,  20),
    "Regional WAN (low)"              : NetworkConfiguration("regWAN100",   100,   20,    15,    0.1,   500,   1200,  10),
    "Cross-continental WAN (high)"    : NetworkConfiguration("ccWAN200",    200,   100,   10,    0.2,   2000,  1200,  20),
    "Cross-continental WAN (medium)"  : NetworkConfiguration("ccWAN100",    100,   120,   15,    0.3,   1000,  1000,  15),
    "Cross-continental WAN (low)"     : NetworkConfiguration("ccWAN50",     50,    150,   25,    0.5,   500,   1000,  10)
}


def set_network_configuration(
    network_configuration,
    port=1212
):
    """
    Configure network settings for local loopback interface with specified parameters.

    This function sets up network emulation on the local loopback interface (`lo`) using 
    `tc` (Traffic Control) commands. The settings include bandwidth limit, delay, jitter, 
    packet loss, burst size, quantum size, and rate to quantum ratio, all of which can
    simulate various network conditions for testing purposes.

    Parameters:
        network_configuration (NetworkConfiguration): a network configuration object.
        port (int, optional): Network port to which these settings apply; default is 1212.

    Raises:
        subprocess.CalledProcessError: If any `tc` command fails to execute.
    """
    print(network_configuration, flush=True)
    commands = [
        # Set up the root qdisc with HTB and a default class
        f"sudo tc qdisc add dev lo root handle 1: htb default 12 {f' r2q {network_configuration.r2q}' if network_configuration.r2q else ''}",
        # Add the class with a rate limit, burst size, quantum, and r2q
        f"sudo tc class add dev lo parent 1: classid 1:1 htb rate {network_configuration.bandwidth}Mbit{f' burst {network_configuration.burst_size}k' if network_configuration.burst_size else ''}{f' quantum {network_configuration.quantum}' if network_configuration.quantum else ''}",
        # Add netem with delay, jitter, and packet loss
        f"sudo tc qdisc add dev lo parent 1:1 handle 10: netem delay {network_configuration.delay}ms{f' {network_configuration.jitter}ms' if network_configuration.jitter else ''}{f' loss {network_configuration.packet_loss}%' if network_configuration.packet_loss else ''}",
        # Apply a filter to limit traffic only to the server port
        f"sudo tc filter add dev lo protocol ip parent 1:0 prio 1 u32 match ip dport {port} 0xffff flowid 1:1"
    ]
    for command in commands:
        subprocess.run(command, shell=True, check=True)


def clean_network_configuration():
    """
    Remove any network emulation settings from the local loopback interface.

    This function deletes all `tc` (Traffic Control) settings on the loopback (`lo`) 
    interface, restoring it to its default state without bandwidth, delay, or packet 
    loss constraints.
    """
    command = "sudo tc qdisc del dev lo root"
    subprocess.run(command, shell=True)


def generate_points(num_dim, num_points, num_clusters):
    cmd = [
        "python3",
        "src/point-generator.py",
        str(num_points),
        str(num_clusters),
        str(num_dim)]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    process.wait()


def get_cmd(party_id, num_dim, num_points, num_clusters, num_iterations):
    filepath = f"data/points_d{num_dim}_n{num_points}_k{num_clusters}.csv"
    if not os.path.exists(filepath):
        print(f"file {filepath} not found -> generating it right now")
        generate_points(num_dim, num_points, num_clusters)
    return["./build/clustering",
           "-k",
           str(num_clusters),
           "-f",
           filepath,
           "-r",
           str(num_iterations),
           "-t",
           "-p",
           str(party_id)]


def run_experiment(num_dim, num_points, num_clusters, num_iterations):

    print(f"Clustering experiment :\n"
          f"    dimensions..{num_dim}\n"
          f"    points......{num_points}\n"
          f"    clusters....{num_clusters}\n"
          f"    rounds......{num_iterations}\n",
          flush=True
         )

    # Define the commands to run
    cmd0 = get_cmd(0, num_dim, num_points, num_clusters, num_iterations)
    cmd1 = get_cmd(1, num_dim, num_points, num_clusters, num_iterations)

    # Start both subprocesses
    process0 = subprocess.Popen(cmd0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    time.sleep(1)
    process1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Collect output and errors from both processes
    output0, error0 = process0.communicate()
    output1, error1 = process1.communicate()

    # Decode output and errors if necessary (e.g., from bytes to string)
    output0 = output0.decode('utf-8').strip()
    output1 = output1.decode('utf-8').strip()
    error0 = error0.decode('utf-8').strip()
    error1 = error1.decode('utf-8').strip()

    # # Check for errors
    # if error0 or error1:
    #     print(error0)
    #     print(error1)
    #     return -1


def repeat_experiment(num_dim, num_points, num_clusters, num_iterations, num_repetitions):
    counter = 0
    while counter < num_repetitions:
        try:
            run_experiment(num_dim, num_points, num_clusters, num_iterations)
            counter += 1
        except:
            pass


def benchmark():
    d           = 2
    t           = 10
    reps        = 10

    clean_network_configuration()

    for network_name, network_configuration in default_network_configurations.items():
        print(network_name)
        set_network_configuration(network_configuration)
        for n in (1000, 10000):
            for k in (2, 5, 8):
                repeat_experiment(d, n, k, t, reps)
        clean_network_configuration()



if __name__ == "__main__":

    # # # # Network configuration parameters
    # # # # network_configuration = default_network_configurations["LAN (low)"]
    # network_configuration = default_network_configurations["LAN (medium)"]
    # # # # network_configuration = default_network_configurations["LAN (high)"]
    # # # # network_configuration = default_network_configurations["Regional WAN (low)"]
    # # # # network_configuration = default_network_configurations["Regional WAN (medium)"]
    # # # # network_configuration = default_network_configurations["Regional WAN (high)"]
    # network_configuration = default_network_configurations["Cross-continental WAN (low)"]
    # # # # network_configuration = default_network_configurations["Cross-continental WAN (medium)"]
    # # # # network_configuration = default_network_configurations["Cross-continental WAN (high)"]
    # # # # network_configuration = default_network_configurations["Corporate WAN (high resilience)"]

    # # Clustering parameters
    # d           = 2
    # n           = 7000
    # k           = 2
    # t           = 1
    # reps        = 1

    # clean_network_configuration()
    # set_network_configuration(network_configuration)
    # repeat_experiment(d, n, k, t, reps)
    # clean_network_configuration()

    benchmark()
