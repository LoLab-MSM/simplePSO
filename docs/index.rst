Welcome to simplepso's documentation!
=====================================

simplePSO is a full implementation of the PSO algorithm. It requires an
objective function to minimize that returns a scalar value, starting position,
bounds on how far you want to search around the starting position, and max
speeds each particle can travel. Then, you set the number of particles in the
swarm, the number of iterations to run, and the number of processors you want
to use. That's it!

.. code-block:: python

    # Here we initial the class
    # We must provide the cost function and a starting value
    optimizer = PSO(start=start_position)

    # We also must set bounds of the parameter space
    optimizer.set_bounds(parameter_range)
    # set speed limits on updating particles
    optimizer.set_speed(speed_min, speed_max)

    # Now we run the pso algorithm
    optimizer.run(num_particles, num_iterations, num_processors,
                  cost_function)
    # get best parameter fit
    best = optimizer.best
    # parameter set
    best.pos
    # fitness of parameter set
    best.fitness

We provide two examples in the examples
`folder <https://github.com/LoLab-VU/simplePSO/tree/master/examples>`_ training
PySB models to time course data using ODE simulations and
`one <https://github.com/LoLab-VU/simplePSO/blob/master/examples/run_schogl_example_ssa.py>`_
example training a model to a distribution of species concentrations using
stochastic simulations.

.. image:: training_earm.gif
   :width: 400

.. toctree::
   :maxdepth: 2

   ode_example/ode_example.rst

.. autoclass:: simplepso.pso.PSO
   :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

