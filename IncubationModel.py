from QuarantineModel import QuarantinedPopulation


class IncubatedPopulation(QuarantinedPopulation):
    def __init__(self, total: int, alpha: float, beta: float, time_step: float, death_rate: float, attendants: int,
                 quarantine_percent: float, symptom_development_rate: float):
        super().__init__(total, alpha, beta, time_step, death_rate, attendants, quarantine_percent)
        self.asymptomatic = self.infected
        self.symptomatic = 0
        self.symptom_development_rate = symptom_development_rate

    def calc_total_attendants(self) -> float:
        """
        :return: total number of susceptible attendants
        """
        return self.susceptible * (1 - (1 - (self.attendants_per_infected / (self.recovered + self.susceptible))) **
                                   (self.quarantine_percent * self.symptomatic))

    def dSdt(self) -> float:
        """
        :return: dS/dt
        """
        return (-self.alpha * 10 * self.symptomatic * self.calc_total_attendants() - self.alpha *
                ((1 - self.quarantine_percent) * self.symptomatic + self.asymptomatic) * self.susceptible) * self.timeStep

    def dIAdt(self) -> float:
        """
        :return: dIA/dt
        """
        return -self.dSdt() - self.symptom_development_rate * self.asymptomatic * self.timeStep

    def dISdt(self) -> float:
        """
        :return: dIS/dt
        """
        return self.symptom_development_rate * self.asymptomatic * self.timeStep - self.beta * self.symptomatic * self.timeStep

    def dRdt(self) -> float:
        """
        :return: dR/dt
        """
        return (1 - self.death_rate) * self.beta * self.symptomatic * self.timeStep

    def dTdt(self) -> float:
        """
        :return: dT/dt
        """
        return self.death_rate * self.beta * self.symptomatic * self.timeStep

    def dSdt_with_masking(self) -> float:
        """
        Highly effective masks for quarantined individuals reduce transmission by 97%
        Cheaper masks reduce general transmission by 50%
        :return: dS/dt
        """
        return (-self.alpha * 0.03 * self.symptomatic * self.calc_total_attendants() - self.alpha * 0.5 *
                ((1 - self.quarantine_percent) * self.symptomatic + self.asymptomatic) * self.susceptible) * self.timeStep

    def dIAdt_with_masking(self) -> float:
        """
        :return: dIA/dt
        """
        return -self.dSdt_with_masking() - self.symptom_development_rate * self.asymptomatic * self.timeStep

    def step(self):
        dSdt = self.dSdt()
        dIAdt = self.dIAdt()
        dISdt = self.dISdt()
        dRdt = self.dRdt()
        dTdt = self.dTdt()

        self.susceptible += dSdt
        self.asymptomatic += dIAdt
        self.symptomatic += dISdt
        self.infected = self.asymptomatic + self.symptomatic
        self.recovered += dRdt
        self.dead += dTdt
        self.time += self.timeStep

    def __str__(self) -> str:
        """
        Gets a string summary of the population demographic
        :return: a summary of the population
        """
        return super(QuarantinedPopulation, self).__str__() + "\nDead: " + str(round(self.dead, 2))
