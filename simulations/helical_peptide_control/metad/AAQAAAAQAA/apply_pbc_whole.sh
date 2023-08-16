for dir in workspace/*; do
    cd $dir;
    for walker_dir in WALKER*; do
        cd $walker_dir;
        gmx trjconv -f npt_new.xtc -s npt_new.tpr -pbc whole -o npt_new.whole.xtc <<< 0;
        cd ..;
    done
    cd ../../;
done

