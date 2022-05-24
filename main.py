from BasicModel import Population
from QuarantineModel import QuarantinedPopulation
from IncubationModel import IncubatedPopulation
import math
import matplotlib.pyplot as plt

# Constants
N = 50000           # Initial Population
alpha = 0.000005    # Transmission Rate
days = 7            # Average Recovery Time After Symptoms
beta = 1/days       # Percentage of Recovery per day
time_step = 1       # Time Step for Simulation
death_rate = 1/10   # Death Rate of Disease
attendants = 0.5      # Number of People that Quarantined Infected Contact
quarantine_percent = 0.7    # Percentage of Population that Quarantines
incubation_period = 3       # Incubation Time Before Symptoms develop
symptomatic_rate = 1/incubation_period      # Percent of Asymptomatic that Become Symptomatic per day
plot_num = 29


def find_quickest_spread(times: [int], susceptible: [float]) -> int:
    """
    Finds when the disease was spreading the fastest
    :param times: a list of times for the simulation
    :param susceptible: susceptible data for the simulation
    :return: the time when the disease was spreading the fastest
    """
    fastest_spread_time = times[0]
    fastest_spread = susceptible[0] - susceptible[1]
    for i in range(1, len(susceptible)-1):
        spread = susceptible[i] - susceptible[i+1]
        if spread > fastest_spread:
            fastest_spread_time = times[i]
            fastest_spread = spread

    return fastest_spread_time


def calc_infected(s: float) -> float:
    """
    Finds the number of infected based on the number of susceptible
    I(S) = beta/alpha * ln(S/InitialS) - S + N
    :param s: the susceptible population
    :return: the infected population
    """
    return beta/alpha * math.log(s/(N-10)) - s + N


def calc_max_infected() -> float:
    """
    Finds the maximum number of people infected at one time
    :return: the maximum number of people infected at one time
    """
    return calc_infected(beta/alpha)


def calc_total_not_infected() -> int:
    """
    Loops through all possible values of the susceptible
    population and finds when I(S) is equal to 0.
    :return: the value of S for which I(S) equals 0.
    """
    for sus in range(N).__reversed__():
        if round(calc_infected(sus)) == 0:
            return sus
    return 0



def main():
    # p = Population(N, alpha, beta, timeStep)
    # p = QuarantinedPopulation(N, alpha, beta, time_step, death_rate, attendants, quarantine_percent)
    p = IncubatedPopulation(N, alpha, beta, time_step, death_rate, attendants, quarantine_percent, symptomatic_rate)
    time_list = [0]
    sus = [p.susceptible]
    inf = [p.infected]
    rec = [p.recovered]
    deaths = []
    if isinstance(p, QuarantinedPopulation):
        deaths.append(p.dead)
    # print(p)
    # print()
    while p.infected > 1:  # Stop simulation once less than half a person is infected
        p.step()  # Simulate one Time Step in Population
        time_list.append(p.time)
        # Add the demographic data of the population
        sus.append(p.susceptible)
        inf.append(p.infected)
        rec.append(p.recovered)
        if isinstance(p, QuarantinedPopulation):
            deaths.append(p.dead)
        # Print out a summary
        print(p)
        print()
    # Find and print when the disease was spreading the fastest
    fast_time = find_quickest_spread(time_list, sus)
    # print(f"{round(rec[-1] + deaths[-1] + inf[-1])}", end=", ")
    print()
    print(f"Total Infected: {round(rec[-1]+deaths[-1]+inf[-1])}")
    print(f"Maximum ill at once: {str(round(max(inf)))} (at time t={str(time_list[inf.index(max(inf))])})")
    print(f"Time of fastest spread: {fast_time} ({str(round(sus[time_list.index(fast_time)] - sus[time_list.index(fast_time) + 1], 2))} p/t)")
    print(f"Total lifespan of disease: {time_list[-1]}")

    # Make a pretty graph
    plt.title("SIR Model " + str(plot_num))
    plt.xlabel("Time")
    plt.ylabel("Number of People")
    plt.plot(time_list, sus, label="Susceptible")
    plt.plot(time_list, inf, label="Infected")
    plt.plot(time_list, rec, label="Recovered")
    if isinstance(p, QuarantinedPopulation) and death_rate > 0.01:
        plt.plot(time_list, deaths, label="Dead")
    plt.legend(loc="upper right")
    plt.axvline(x=fast_time, color="red")
    # plt.text(fast_time - 8.5, 0, "t=" + str(fast_time))
    plt.savefig("SIR_Model_" + str(plot_num))
    plt.show()


if __name__ == "__main__":
    # plt.title("Altering Incubation Period")
    # plt.xlabel("Incubation Period")
    # plt.ylabel("Total Number of Infected")
    # nums = [42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 49, 49, 50, 50, 51, 52, 53, 53, 54, 55, 56, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 69, 70, 71, 73, 74, 75, 77, 79, 80, 82, 84, 86, 88, 90, 92, 94, 96, 99, 102, 104, 107, 110, 114, 117, 121, 125, 129, 133, 138, 143, 149, 155, 161, 168, 176, 184, 193, 203, 214, 227, 240, 256, 273, 293, 315, 341, 370, 405, 445, 492, 548, 614, 693, 785, 893, 1018, 1160, 1318, 1491, 1679, 1879, 2089, 2309, 2537, 2771, 3012, 3259, 3510, 3765, 4025, 4288, 4555, 4825, 5098, 5375, 5654, 5936, 6221, 6509, 6799, 7092, 7387, 7685, 7985, 8287, 8592, 8899, 9208, 9519, 9832, 10147, 10464, 10783, 11103, 11425, 11749, 12074, 12400, 12727, 13056, 13385, 13716, 14047, 14379, 14711, 15043, 15376, 15708, 16041, 16373, 16706, 17037, 17368, 17698, 18027, 18355, 18681, 19007, 19330, 19653, 19973, 20292, 20608, 20922, 21234, 21544, 21852, 22156, 22459, 22758, 23055, 23349, 23640, 23927, 24212, 24494, 24773, 25049, 25321, 25590, 25856, 26119, 26378, 26634, 26887, 27136, 27382, 27625, 27865, 28101, 28334, 28564, 28791, 29014, 29235, 29452, 29666, 29877, 30085, 30290, 30492, 30691, 30887, 31080, 31270, 31458, 31642, 31824, 32004, 32180, 32354, 32526, 32695, 32861, 33025, 33187, 33346, 33503, 33658, 33810, 33960, 34108, 34254, 34398, 34539, 34679, 34817, 34952, 35086, 35218, 35348, 35476, 35602, 35727, 35849, 35970, 36090, 36208, 36324, 36438, 36551, 36663, 36773, 36881, 36988, 37094, 37198, 37301, 37403, 37503, 37602, 37699, 37796, 37891, 37985, 38078, 38169, 38259, 38349, 38437, 38524, 38610, 38695, 38779, 38862, 38944, 39024, 39104, 39183, 39261, 39338, 39415, 39490, 39564, 39638, 39711, 39783, 39854, 39924, 39993, 40062, 40130, 40197, 40264, 40329, 40394, 40459, 40522, 40585, 40647, 40709, 40770, 40830, 40889, 40948, 41007, 41065, 41122, 41178, 41234, 41290, 41345, 41399, 41453, 41506, 41558, 41611, 41662, 41713, 41764, 41814, 41864, 41913, 41962, 42010, 42058, 42105, 42152, 42198, 42244, 42290, 42335, 42380, 42424, 42468, 42511, 42554, 42597, 42639, 42681, 42723, 42764, 42805, 42845, 42886, 42925, 42965, 43004, 43043, 43081, 43119, 43157, 43194, 43231, 43268, 43305, 43341, 43377, 43412, 43448, 43483, 43517, 43552, 43586, 43620, 43653, 43687, 43720, 43752, 43785, 43817, 43849, 43881, 43913, 43944, 43975, 44006, 44036, 44067, 44097, 44126, 44156, 44185, 44215, 44244, 44272, 44301, 44329, 44357, 44385, 44413, 44440, 44468, 44495, 44522, 44548, 44575, 44601, 44627, 44653, 44679, 44704, 44730, 44755, 44780, 44805, 44829, 44854, 44878, 44902, 44926, 44950, 44974, 44997, 45021, 45044, 45067, 45090, 45112, 45135, 45157, 45180, 45202, 45224, 45246, 45267, 45289, 45310, 45331, 45352, 45373, 45394, 45415, 45436, 45456, 45476, 45496, 45516, 45536, 45556, 45576, 45595, 45615, 45634, 45653, 45672, 45691, 45710, 45729, 45747, 45766, 45784, 45802, 45821, 45839, 45857, 45874, 45892, 45910, 45927, 45945, 45962, 45979, 45996, 46013, 46030, 46047, 46063, 46080, 46096, 46113, 46129, 46145, 46161, 46177, 46193, 46209, 46225, 46240, 46256, 46272, 46287, 46302, 46317, 46332, 46348, 46362, 46377, 46392, 46407, 46421, 46436, 46450, 46465, 46479, 46493, 46507, 46522, 46535, 46549, 46563, 46577, 46591, 46604, 46618, 46631, 46645, 46658, 46671, 46684, 46698, 46711, 46724, 46736, 46749]
    # hi = list(range(100, 600))
    # for i in range(len(hi)):
    #     hi[i] /= 100
    # # plt.plot(hi, nums, "r*")
    # plt.plot(hi, nums)
    # #
    # plt.savefig("SIR_Model_21")
    # plt.show()
    for i in range(1):
        main()




    # print(calc_total_not_infected())