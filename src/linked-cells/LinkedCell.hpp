#pragma once

#include <list>
#include <memory>

#include "../utility/utility.hpp"

using Particle = Utility::Particle;
class LinkedCell {
private:
    std::list<std::shared_ptr<Particle>> particles_;

public:
    LinkedCell() {}

    void addParticle(const std::shared_ptr<Particle>& particle) { particles_.push_back(particle); }

    void removeParticle(const std::shared_ptr<Particle>& particle) { particles_.remove(particle); }

    [[nodiscard]] const std::list<std::shared_ptr<Particle>>& getParticles() const { return particles_; }
};