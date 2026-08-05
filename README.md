# UNISIS

### Clone

```
git clone git://github.com/naotohori/sis.git
```

Required submodules will be automatically installed via cmake when build. If it does not work, use the following git command.

```
git clone --recurse-submodules -j8 git://github.com/naotohori/sis.git
```

### Compile

```
mkdir -p ./build && cd ./build
cmake ..
make -j12
cd ..
ln -s build/sis
```
