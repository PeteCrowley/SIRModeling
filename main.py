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
attendants = 5      # Number of People that Quarantined Infected Contact
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
    main()
