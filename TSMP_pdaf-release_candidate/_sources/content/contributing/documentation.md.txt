# Contribute to TSMP documentation.

Contributing to the TSMP documentation should be achieve via a pull request from 
a personal fork of the TSMP repo. That way no unreviewed changes will find their 
way into the official documentation, makes it easy to discuss upcoming changes 
via related issues, and provides a clear to-do list for all open pull requests. 

In detail you need to fork the TSMP repo, create a local clone of the fork, and 
update / change the documentation locally. When you are finished, commit your 
changes with a meaningful commit message, push the new commits to your remote 
fork and create a pull request. 

For your local changes, all you need is a text editor to modify the individual 
`.md` files that make up the documentation. However, to see if your changes 
render as expected, we recommend that you use sphinx to build the entire 
documentation. Do not use a markdown editor for this, as the markdown flavour 
used may differ from that used by sphinx, leading to different results. 

The internal structure of the documentation is quite simple. All files related 
to the documentation are located in `doc/`.   
`conf.py` controls the behaviour of sphinx, `index.rst` is the entry point to 
the documentation, and `content/` contains all the individual `.md`files.  
Building the documentation is done by moving to `doc/` and running the command: 
```
sphinx-build -a . _build
``` 
This will create the documentation in `doc/_build/`. 
Simply browse to this directory and open `index.html` which should show you the 
locally rendered documentation in your default web browser.

