class TrackerBeamline:

    def __init__(self, trackers = []):
        self.trackers = trackers
        self.particles = []
        self.input_particles = []

    def set_input_particles(self, particles):
        self.input_particles = particles

    def add_tracker(self, tracker):
        self.trackers.append(tracker)

    def track(self, save_step = False ):
        # empty particles
        self.particles = []

        # initial particle distribution
        particles = self.input_particles

        # loop and apply trackers in order
        for t in self.trackers:

            if save_step:
                self.particles.append(particles)

            # track through tracker and update particles
            particles = t.track(particles)


        return particles
