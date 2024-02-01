#pragma once

#include <array>
#include <string>


class Particle {

private:


    /**
     * Force effective on this particle
     */
    std::array<double, 3> f{};

    /**
     * Force which was effective on this particle
     */
    std::array<double, 3> old_f{};

    /**
     * Position of the particle
     */
    std::array<double, 3> x{};

    /**
     * Velocity of the particle
     */
    std::array<double, 3> v{};

    /**
     * Mass of this particle
     */
    double m{};

    std::array<double, 3> grav{};

    /**
     * the base force applied in each iteration
     */
    std::array<double, 3> baseForce{};

    /**
     * original index of the particle in the grid generated. used for membranes
     */
    std::array<int, 3> gridIndex{};

    std::array<double, 3> oldX{};

    bool passedPeriodicX = false;
    bool passedPeriodicY = false;
    bool passedPeriodicZ = false;

public:
    [[nodiscard]] bool isPassedPeriodicY() const;

    [[nodiscard]] bool isPassedPeriodicX() const;

    [[nodiscard]] bool isPassedPeriodicZ() const;

    /**
     * Whether this particle is still valid and should be printed/used for calculations
     */
    bool valid;


    /**
     * whether the Particle is part of a membrane -> different calculation
     */
    bool membrane;

    /**
     * @brief Default constructor
     */
    explicit Particle();


    /**
     * @brief Copy constructor
     * @param other object it is copied from
     */
    Particle(const Particle &other);

    /**
     * @brief Copy assignment operator
     * @param other object it is copied from
     * @return returning the copy
     */
    Particle &operator=(const Particle &other);

    /**
     * @brief Move constructor
     * @param other object it is copied from
     */
    Particle(Particle &&other) noexcept;

    /**
     * @brief Move assignment operator
     * @param other object it is copied from
     * @return returning the copy
     */
    Particle &operator=(Particle &&other) noexcept;

    bool operator<(const Particle &other) noexcept;

    /**
     * @brief Overloaded constructor
     * @param x_arg position
     * @param v_arg velocity
     * @param m_arg mass
     */
    Particle(
            // for visualization, we need always 3 coordinates
            // -> in case of 2d, we use only the first and the second
            std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg);


    virtual ~Particle();

    [[nodiscard]] const std::array<double, 3> &getX() const;

    [[nodiscard]] const std::array<double, 3> &getV() const;

    [[nodiscard]] const std::array<double, 3> &getF() const;

    [[nodiscard]] const std::array<double, 3> &getOldF() const;

    [[nodiscard]] double getM() const;

    [[nodiscard]] const std::array<double, 3> &getGrav() const;


    [[nodiscard]] const std::array<int, 3> &getGridIndex() const;

    [[nodiscard]] const std::array<double, 3> &getOldX() const;

    [[nodiscard]] const std::array<double, 3> &getBaseForce() const;


    void setPassedPeriodic(int d);

    void setMembrane(bool membrane);

    void setGridIndex(std::array<int, 3> &index);


    /**
     * sets the base force for a particle
     * @param force the base force
     */
    void setBaseForce(std::array<double, 3> &force);

    /**
     * applies the base force and a given gravity.
     * sets force of particle to baseForce attribute and add grav*mass in y direction
     * @param grav the gravity factor
     */
    void applyBaseForceAndGrav();

    /**
     * @param d: dimension
     * @param v: value
     */
    void setX(int d, double v);

    /**
     * set the position vector to some other vector
     * @param x the new position vector
     */
    void setX(std::array<double, 3> x);

    void setGrav(std::array<double, 3> grav);

    void setOldX(std::array<double, 3> x);

    /**
     * @param d: dimension
     * @param v: value
     */
    void setV(int d, double v);

    /**
     * set the velocity vector to some other vector
     * @param v the new velocity vector
     */
    void setV(std::array<double, 3> v);

    /**
     * @param d: dimension
     * @param v: value
     */
    void setF(int d, double v);

    /**
     * set the force vector to some other vector
     * @param f the new force vector
     */
    void setF(std::array<double, 3> f);

    /**
     * @param val the old force array
     */
    void setOldF(const std::array<double, 3> &val);

    void setM(double);


    /**
     * @brief Class-specific equals operator
     * @param other object it is compared to
     * @return true if equal, false if not
     */
    bool operator==(Particle &other);

    /**
     * @brief Class-specific not-equals operator
     * @param other object it is compared to
     * @return true if not equal, false if equal
     */
    bool operator!=(Particle &other);

    /**
     * @brief Turns the particle into a string
     * @return the string
     */
    [[nodiscard]] std::string toString() const;


};


//void to_json(json&, const Particle&);

//void from_json(const json&, Particle& p);


/**
 * @brief Overloaded << operator for easy console prints
 * @param stream the stream the content is pushed into into
 * @param p the particle
 * @return the stream reference for further pushing into to stream
 */
std::ostream &operator<<(std::ostream &stream, Particle &p);
