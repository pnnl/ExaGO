FROM ubuntu:18.04

RUN apt update && apt install -y python3 python3-pip \
  gpg git gcc-7 g++-7 gfortran-7 file findutils automake autoconf make \
  patch tar xz-utils gzip unzip bzip2 zstd bash

RUN touch /cache.exasgd.pnl.gov.pub \
&& (echo "-----BEGIN PGP PUBLIC KEY BLOCK-----" \
&& echo "" \
&& echo "mQINBGEwBrABEADQ+aBOACIXLJJaJnZK5XjBCMcIyPDqZ6Lsw0Es88gMXEb+RZRQ" \
&& echo "hs+UyLoc4hFn/CnvXM4nqCv+nz2cMyWz6zttrrc9JMuF1PE16Fg0qmiA7I8tTLBK" \
&& echo "lg99lKQwypSP7Mers3IX7MT13WwTXn6J7wPeguke9O8ZQqJlC4Q8TRG7KtmkRDpj" \
&& echo "FSC8t+ah0fdQ9UB3oiJJAes1i4bNkIfrPvtswUfO8v7s+Gp4Elg1T9P+W/XwII7k" \
&& echo "3EReN9+AZegol2FNUfFsIYN8HWX4qk80uLXdJ6p/hxw2PKTq6YMKLClFtpHVzEY8" \
&& echo "YC6Y/2bGqfqg66kdsOyRGcaxSu5gSKWT4KnnFCnVNPUvz0FVdRzwcbuOCRxiLk+7" \
&& echo "kPA2UnXOH1mcraDp1zTiuILKHlybtly1oOslgIvIwUG04iJm5wjrVKn5n0NV1nn3" \
&& echo "Wk78slFqdmBnTh5rEjBH6IM+CshYWz4c54YhG0kq5WpUfGvkyHFh+gh6jKfHexLN" \
&& echo "JMWLgm+pWo58pTNtuQwIxA3eEr0W+leC9Un21oZBKvwhZXqu/pLmJI2FBUxq/IE6" \
&& echo "ZuxUZKuJJmzbt54CIzwZnldxrd2hFiN0Zbgdm75MCbvk7xbcGAoT5byf0MXhHquN" \
&& echo "Pp/l3/rYYHI9vEDrPNxuZFEWsaOeZuwasH15VwlFw1OAFQWiNNJcnKU7XQARAQAB" \
&& echo "tERBc2hlciBNYW5jaW5lbGxpIChHUEcgY3JlYXRlZCBmb3IgU3BhY2spIDxhc2hl" \
&& echo "cm1hbmNpbmVsbGlAZ21haWwuY29tPokCTgQTAQoAOBYhBPmj4IaICYKDdFUUEuHj" \
&& echo "WhxvhNadBQJhMAawAhsDBQsJCAcCBhUKCQgLAgQWAgMBAh4BAheAAAoJEOHjWhxv" \
&& echo "hNadiMQP/jI8/ATsAVP0xUW287dbGGI2ypsllJzRR800woXreIEJEsXyZr+GR8La" \
&& echo "N7/jjKcfgBkTbSGaX1jNgtG/gcPojDGe+KqbHTItQKxjhxD1Bo1gQUJ5UPqYkSFc" \
&& echo "BruF9fMIJnA89b97gqBMLc6Of1uq53EP8VQDCq1Vg0J8j7sbL7CYldk6rila8gZF" \
&& echo "dFK8DBkUMw79MRE9YizIGfu6XYNtTf5TXrlA1SRRgEtkyjbQ1i7qsgu7e9MR+yi7" \
&& echo "yD1Rp3HxgPmgA9AEXGWCSAvooo9phzocXp1tVVZa+LU5fT73Wras182MzdUzQoyw" \
&& echo "SYRWnBqIZajLymwdUZoNfOWQ65nWgMFRwP0wCgs8e4GknYM46Bz84DUXmTPCkjow" \
&& echo "aGZ9xIEi5apsdpzuEPsH1tMr0A/0qwGX3qw4xeY5Q+nLnn2/5gEnSauLbjQjAkIt" \
&& echo "LfEPa36k57cCszyleyx+gtoACx5fvgtUUEoQdJFIiLyLY1uIbE9ZjOlGRb5cC/Xe" \
&& echo "Bo1fQA7TpuuzXeelxa4CK4GiidGysGCLdZxf9vZqC4um8NscYDcL+N+XJtrs5pht" \
&& echo "QYAmjWfDdfN1fg/Bl2w4tpUqVhd4Wva5eGJFfxo++6ijnmL3TPh3BBdGLzc1Kcfr" \
&& echo "3FaOTIYCoopuhzjU+xLnhRp0ynuLI6ZTTKjc+TjRQeIgrbevfPQh" \
&& echo "=1iQq" \
&& echo "-----END PGP PUBLIC KEY BLOCK-----") > /cache.exasgd.pnl.gov.pub

RUN mkdir -p /opt/spack-environment \
&&  (echo "spack:" \
&&   echo "  specs:" \
&&   echo "  - exago@develop~ipopt+hiop+mpi~cuda ^hiop@0.4.1 ^openmpi ^petsc@3.13.5+mpi~hypre~hdf5~metis~superlu-dist" \
&&   echo "  config:" \
&&   echo "    clingo: true" \
&&   echo "    install_tree: /opt/software" \
&&   echo "  concretization: together" \
&&   echo "  view: /opt/view" \
&&   echo "  mirrors:" \
&&   echo "    minio: s3://spack") > /opt/spack-environment/spack.yaml

RUN git clone -b develop --single-branch https://github.com/LLNL/spack.git /spack

RUN pip3 install boto boto3 botocore

# Install the software, remove unnecessary deps
ENV S3_ENDPOINT_URL=http://cache.exasgd.pnl.gov
RUN cd /opt/spack-environment && \
  . /spack/share/spack/setup-env.sh && \
  spack env activate . && \
  spack gpg init && \
  spack gpg trust /cache.exasgd.pnl.gov.pub && \
  spack buildcache update-index -d s3://spack && \
  spack install --fail-fast

