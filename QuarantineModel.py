from BasicModel import Population


class QuarantinedPopulation(Population):
    def __init__(self, total: int, alpha: float, beta: float, time_step: float, death_rate: float, attendants: int,
                 quarantine_percent: float):
        super().__init__(total, alpha, beta, time_step)
        self.dead = 0
        self.death_rate = death_rate
        self.attendants_per_infected = attendants
        self.quarantine_percent = quarantine_percent

    def calc_total_attendants(self) -> float:
        return self.susceptible*(1-(1-(self.attendants_per_infected/(self.recovered + self.susceptible))) **
                                    (self.quarantine_percent * self.infected))

    def dSdt(self) -> float:
        """
        :return: dS/dt
        """
        return (-self.alpha * self.infected * self.calc_total_attendants() - self.alpha *
                (1-self.quarantine_percent) * self.infected * self.susceptible) * self.timeStep

    def dIdt(self) -> float:
        """
        :return: dI/dt
        """
        return -self.dSdt() - self.beta * self.infected * self.timeStep

    def dRdt(self) -> float:
        """
        :return: dR/dt
        """
        return (1-self.death_rate)*self.beta * self.infected * self.timeStep

    def dTdt(self) -> float:
        """
        :return: dT/dt
        """
        return self.death_rate * self.beta * self.infected * self.timeStep

    def step(self):
        dSdt = self.dSdt()
        dIdt = self.dIdt()
        dRdt = self.dRdt()
        dTdt = self.dTdt()

        self.susceptible += dSdt
        self.infected += dIdt
        self.recovered += dRdt
        self.dead += dTdt
        self.time += self.timeStep

    def __str__(self) -> str:
        """
        Gets a string summary of the population demographic
        :return: a summary of the population
        """
        return super(QuarantinedPopulation, self).__str__() + "\nDead: " + str(round(self.dead, 2))








