// ============================================================================
//  fielddump -- flatten one variable of a single-level 2-D AMReX plotfile onto
//  a regular array, for the NSCBC multi-dimensional metrics.
//
//  (Named fielddump rather than pltdump: the repository .gitignore carries a
//  `plt*` rule for plotfiles, which silently swallows any source file whose
//  name begins with "plt".)
//
//  The 2-D boundary tests need the whole field, not a line-out, because the
//  quantity of interest is the SHAPE of the wavefront.  fextract gives 1-D
//  slices only, so this exists.
//
//    fielddump <plotfile> <varname> <outfile>
//
//  Writes an ASCII header (nx ny xlo ylo dx dy time) followed by nx*ny values
//  in row-major (j slowest) order.  Deliberately trivial to parse.
// ============================================================================

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>

#include <cstdio>
#include <string>
#include <vector>

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, false);
  int rc = 0;
  {
    if (argc < 4) {
      amrex::Print() << "usage: fielddump <plotfile> <varname> <outfile>\n";
      amrex::Finalize();
      return 1;
    }
    const std::string pf(argv[1]);
    const std::string var(argv[2]);
    const std::string out(argv[3]);

    amrex::PlotFileData plotfile(pf);
    const int lev = 0; // these tests are single level by construction
    const amrex::Box dom = plotfile.probDomain(lev);
    const auto plo = plotfile.probLo();
    const auto dx = plotfile.cellSize(lev);

    bool found = false;
    for (const auto& n : plotfile.varNames()) {
      found = found || (n == var);
    }
    if (!found) {
      amrex::Print() << "fielddump: variable '" << var << "' not in " << pf
                     << "\n  available:";
      for (const auto& n : plotfile.varNames()) {
        amrex::Print() << " " << n;
      }
      amrex::Print() << "\n";
      amrex::Finalize();
      return 1;
    }

    const amrex::MultiFab mf = plotfile.get(lev, var);

    const int nx = dom.length(0);
    const int ny = dom.length(1);
    std::vector<double> a(
      static_cast<size_t>(nx) * ny, std::numeric_limits<double>::quiet_NaN());

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.validbox();
      const auto& arr = mf.const_array(mfi);
      const auto lo = amrex::lbound(bx);
      const auto hi = amrex::ubound(bx);
      for (int j = lo.y; j <= hi.y; ++j) {
        for (int i = lo.x; i <= hi.x; ++i) {
          a[static_cast<size_t>(j - dom.smallEnd(1)) * nx +
            (i - dom.smallEnd(0))] = arr(i, j, lo.z);
        }
      }
    }

    FILE* f = std::fopen(out.c_str(), "w");
    std::fprintf(
      f, "# nx ny xlo ylo dx dy time\n%d %d %.17g %.17g %.17g %.17g %.17g\n",
      nx, ny, plo[0], plo[1], dx[0], dx[1], plotfile.time());
    for (int j = 0; j < ny; ++j) {
      for (int i = 0; i < nx; ++i) {
        std::fprintf(f, "%.10e ", a[static_cast<size_t>(j) * nx + i]);
      }
      std::fprintf(f, "\n");
    }
    std::fclose(f);
    amrex::Print() << "fielddump: wrote " << nx << " x " << ny << " '" << var
                   << "' at t = " << plotfile.time() << " to " << out << "\n";
  }
  amrex::Finalize();
  return rc;
}
