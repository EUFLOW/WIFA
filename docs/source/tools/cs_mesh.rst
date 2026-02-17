Geometry and Mesh
--------------------
The model setup embedded within WIFA includes a circular computational domain centered on the farm automatically generated using the open-source SALOME platform (python scripts included in WIFA). The default domain is three times bigger than the farm and 800 m high. The minimal mesh element size is equal to D/12 at turbine positions and cell size progressively becomes larger away from the farm up to 0.8×D at lateral boundaries. These characteristics can be modified by the user through the python API.
