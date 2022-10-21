---
title: Helpnote 1 - Jupyter on HEP Cluster
author: Till Dieminger
---

We use the following terms:

- Server Side: The Manchester HEP cluster at `<hostname>`
- Client Side: Your local machine

## Set up python

For this we use the [`cvmfs` library](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html) which is already installed and mounted on the HEP cluster.
To activate what we need, log on to the server side and run

```sh
source /cvmfs/larsoft.opensciencegrid.org/products/setup

setup python v3_9_2
```
As we will need serveral libraries, now is a good time to check if this worked and do some setup at the same time.
Run

```sh
pip install numpy scipy matplotlib ipython jupyter pandas uproot awkward
```
This should now work - if not, lets see.

## Start the notebook

On the server side execute
```sh
python -m notebook --no-browser --port=7800
```
This should fire up the notebook on the server side.
In the upcoming dialoge there will be a note giving you the token you need later.



On the client side now run

```sh
ssh -N -f -L localhost:8001:localhost:7800 <username>@<hostname>
```

where `<username>` and `<hostname>` should be replaced by your corresponding credentials and the hostname of the HEP cluster you work on.

Now open your browser and navigate to `http://localhost:8001/`.
Here the notebook should appear.
Enter the token mentioned before from the server side and you should be up and running.
