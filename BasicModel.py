import matplotlib.pyplot as plt

class Population:
    def __init__(self, total: int, alpha: float, beta: float, timeStep: float):
        self.infected = 10                          # Initial Infected Population
        self.susceptible = total - self.infected    # Initial Susceptible Population
        self.recovered = 0                          # Initial Recovered Population
        self.time = 0                               # Initial Time
        self.alpha = alpha
        self.beta = beta
        self.timeStep = timeStep

    def dSdt(self) -> float:
        """
        dS/dt = -alpha * Infected * Susceptible
        :return: dS/dt
        """
        return -self.alpha * self.infected * self.susceptible * self.timeStep

    def dIdt(self) -> float:
        """
        dI/dt = alpha * Infected * Susceptible - beta * Infected
        :return: dI/dt
        """
        return -self.dSdt() - self.dRdt()

    def dRdt(self) -> float:
        """
        dR/dt = beta * Infected
        :return: dR/dt
        """
        return self.beta * self.infected * self.timeStep

    def step(self):
        dSdt = self.dSdt()
        dIdt = self.dIdt()
        dRdt = self.dRdt()

        self.susceptible += dSdt
        self.infected += dIdt
        self.recovered += dRdt

        self.time += self.timeStep

    def __str__(self) -> str:
        """
        Gets a string summary of the population demographic
        :return: a summary of the population
        """
        return "Time: " + str(self.time) + "\nSusceptible: " + str(round(self.susceptible, 2)) + "\nInfected: " + \
               str(round(self.infected, 2)) + "\nRecovered: " + str(round(self.recovered, 2))










